#!/usr/bin/env python3
"""
main.py — Physics-Informed 3D Graph ML для прогноза каталитической посадки субстрата.

Особенности:
  • Многоуровневое логирование (консоль + outputs/logs/main.log).
  • Валидация PDB до вызова Bio.PDB.PDBParser (отсеивает HTML/gzip/mmCIF).
  • Явные исключения с осмысленными сообщениями.
  • Graceful fallback на 3D-граф субстрата, если структура фермента битая.
  • Замер и вывод таймингов + окружения (важно для WSL).
"""

from __future__ import annotations

import sys
import time
import platform
import logging
import traceback
from pathlib import Path

from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR

# =========================================================
# 1. Импорт тяжёлых зависимостей с явной диагностикой
# =========================================================
_DEP_ERRORS: list[str] = []

try:
    import pandas as pd
except ImportError as e:
    _DEP_ERRORS.append(f"pandas: {e}")

try:
    import torch
    import torch.nn as nn
except ImportError as e:
    _DEP_ERRORS.append(f"torch: {e}")

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, AllChem

    RDLogger.DisableLog("rdApp.*")  # убираем шум RDKit в stdout
except ImportError as e:
    _DEP_ERRORS.append(f"rdkit: {e}")

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionException
except ImportError as e:
    _DEP_ERRORS.append(f"biopython: {e}")

if _DEP_ERRORS:
    sys.stderr.write(
        "\n❌ Не хватает Python-зависимостей:\n  - "
        + "\n  - ".join(_DEP_ERRORS)
        + "\n\nУстанови их и повтори запуск. Пример для conda:\n"
        "  pip install torch pandas rdkit biopython\n\n"
    )
    sys.exit(2)


# =========================================================
# 2. Логирование (консоль + файл)
# =========================================================
LOG_DIR = OUTPUT_DIR / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOG_DIR / "main.log"

_fmt = "%(asctime)s [%(levelname)-7s] %(name)-18s | %(message)s"
_datefmt = "%H:%M:%S"

logging.basicConfig(
    level=logging.INFO,
    format=_fmt,
    datefmt=_datefmt,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, mode="a", encoding="utf-8"),
    ],
)
logger = logging.getLogger("Physics3DGraph")


# =========================================================
# 3. Константы путей и кастомные исключения
# =========================================================
DATASET_PATH = DATA_DIR / "CatPred-DB/data/kcat/kcat-random_trainvaltest.csv"
PDB_DIR = DATA_DIR / "alphafold_pdb"


class PipelineError(Exception):
    """Общий предок ошибок пайплайна."""


class DatasetError(PipelineError):
    """Проблема с CSV CatPred-DB."""


class LigandError(PipelineError):
    """SMILES не парсится / RDKit не смог построить 3D."""


class ProteinFormatError(PipelineError):
    """Файл есть, но это не PDB (HTML, gzip, mmCIF, пустой)."""


class ProteinParseError(PipelineError):
    """PDBParser упал или не нашёл ни одного CA."""


# =========================================================
# 4. Модель
# =========================================================
class PhysicsAdd3DGraphNet(nn.Module):
    """Эквивариантная физически-дополненная сеть оценки энергии посадки."""

    def __init__(self):
        super().__init__()
        self.parser = PDBParser(QUIET=True)
        self.graph_energy_head = nn.Sequential(
            nn.Linear(4, 32),
            nn.Tanh(),
            nn.Linear(32, 1),
        )
        logger.debug("[model] инициализирован PhysicsAdd3DGraphNet")

    # ----------------------------------------------------
    def _validate_pdb_format(self, path: Path) -> None:
        """
        Проверяет, что path — действительно PDB. Бросает ProteinFormatError
        с точным диагнозом, если файл HTML/gzip/mmCIF/пустой.
        """
        if not path.exists():
            raise ProteinFormatError(f"файл не существует: {path}")
        size = path.stat().st_size
        logger.info("[validate] размер файла: %d байт", size)
        if size == 0:
            raise ProteinFormatError(f"файл пустой: {path}")

        with open(path, "rb") as f:
            head = f.read(128)
        stripped = head.lstrip().lower()

        if stripped.startswith(b"<!doctype html") or stripped.startswith(b"<html"):
            raise ProteinFormatError(
                f"{path} — HTML-страница, а не PDB. "
                f"Скачай прямую ссылку AF-<UNIPROT>-F1-model_v4.pdb "
                f"(scripts/fetch_alphafold.sh) или проверь VPN/DNS."
            )
        if head[:2] == b"\x1f\x8b":
            raise ProteinFormatError(
                f"{path} — gzip. Распакуй: gunzip -c {path}.gz > {path}"
            )
        if head.startswith(b"data_"):
            raise ProteinFormatError(
                f"{path} — mmCIF (начинается с 'data_'). Нужен MMCIFParser "
                f"или перекачай PDB-формат."
            )

        n_atom = 0
        with open(path, "r", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    n_atom += 1
        if n_atom == 0:
            raise ProteinFormatError(f"{path}: 0 строк ATOM/HETATM — это не PDB.")
        logger.info("[validate] формат OK: %d ATOM/HETATM строк", n_atom)

    # ----------------------------------------------------
    def extract_protein_coords(self, pdb_path: Path) -> torch.Tensor:
        """
        Возвращает [N, 3] float32 тензор координат CA-атомов.
        Бросает ProteinFormatError / ProteinParseError.
        """
        logger.info("[extract] старт: %s", pdb_path)
        self._validate_pdb_format(pdb_path)

        try:
            structure = self.parser.get_structure("protein", str(pdb_path))
        except (PDBConstructionException, ValueError) as e:
            raise ProteinParseError(f"PDBParser упал на {pdb_path}: {e}") from e

        coords = []
        n_models = n_chains = n_res = n_ca = 0
        for model in structure:
            n_models += 1
            for chain in model:
                n_chains += 1
                for residue in chain:
                    n_res += 1
                    if "CA" in residue:
                        n_ca += 1
                        coords.append(residue["CA"].get_coord())

        logger.info(
            "[extract] models=%d chains=%d residues=%d CA=%d",
            n_models,
            n_chains,
            n_res,
            n_ca,
        )
        if not coords:
            raise ProteinParseError(
                f"в {pdb_path} не найдено ни одного CA "
                f"(residues={n_res}). Проверь, что это полноценный PDB."
            )

        tensor = torch.tensor(coords, dtype=torch.float32)
        logger.info(
            "[extract] ✓ shape=%s ndim=%d dtype=%s",
            tuple(tensor.shape),
            tensor.ndim,
            tensor.dtype,
        )
        logger.debug("[extract] первые 3 CA: %s", tensor[:3].tolist())
        return tensor

    # ----------------------------------------------------
    def forward(self, protein_coords: torch.Tensor, tpsa: float, logp: float):
        logger.debug(
            "[forward] in: shape=%s ndim=%d numel=%d dtype=%s",
            tuple(protein_coords.shape),
            protein_coords.ndim,
            protein_coords.numel(),
            protein_coords.dtype,
        )

        if protein_coords.numel() == 0:
            raise ValueError("protein_coords пуст")
        if protein_coords.ndim == 1 and protein_coords.shape[0] == 3:
            logger.warning("[forward] protein_coords 1D [3] → unsqueeze → [1, 3]")
            protein_coords = protein_coords.unsqueeze(0)
        if protein_coords.ndim != 2 or protein_coords.shape[1] != 3:
            raise ValueError(
                f"protein_coords неверной формы {tuple(protein_coords.shape)}, ожидалось [N, 3]"
            )

        w_barrier = (tpsa * 0.4) - (logp * 1.2)
        logger.info("[forward] w_barrier = %.4f кДж/моль", w_barrier)

        protein_center = torch.mean(protein_coords, dim=0)
        logger.info(
            "[forward] center = [%.2f, %.2f, %.2f]",
            *protein_center.tolist(),
        )

        distances = torch.linalg.vector_norm(
            protein_coords - protein_center, dim=-1, keepdim=True
        )
        n = int(protein_coords.shape[0])
        logger.info(
            "[forward] distances: shape=%s (n_nodes=%d)", tuple(distances.shape), n
        )

        tpsa_feat = torch.full((n, 1), float(tpsa), dtype=torch.float32)
        logp_feat = torch.full((n, 1), float(logp), dtype=torch.float32)
        barrier_feat = torch.full((n, 1), float(w_barrier), dtype=torch.float32)

        node_features = torch.cat(
            [distances, tpsa_feat, logp_feat, barrier_feat], dim=1
        )
        logger.info("[forward] node_features: shape=%s", tuple(node_features.shape))

        node_energies = self.graph_energy_head(node_features)
        pose_prob = torch.softmax(-node_energies, dim=0)
        effective_aff = torch.sum(pose_prob * node_energies)
        pose_probability = torch.sigmoid(effective_aff)

        logger.info(
            "[forward] ✓ affinity=%.4f pose_prob=%.4f",
            effective_aff.item(),
            pose_probability.item(),
        )
        return pose_probability, protein_center, w_barrier


# =========================================================
# 5. Вспомогательные функции пайплайна
# =========================================================
def log_environment() -> None:
    logger.info("=" * 60)
    logger.info("Окружение:")
    logger.info("  python  : %s", sys.version.split()[0])
    logger.info("  platform: %s", platform.platform())
    logger.info("  cwd     : %s", Path.cwd())
    logger.info("  torch   : %s", torch.__version__)
    logger.info("  cuda    : %s", "да" if torch.cuda.is_available() else "нет")
    is_wsl = "microsoft" in platform.uname().release.lower()
    logger.info("  WSL     : %s", "да" if is_wsl else "нет")
    logger.info("  log file: %s", LOG_FILE.resolve())
    logger.info("=" * 60)


def load_first_row(csv_path: Path) -> tuple[str, str]:
    """Читает первую строку датасета. Возвращает (uniprot, smiles)."""
    if not csv_path.exists():
        raise DatasetError(f"датасет не найден: {csv_path.resolve()}")
    try:
        df = pd.read_csv(csv_path, nrows=1)
    except Exception as e:
        raise DatasetError(f"pandas не смог прочитать {csv_path}: {e}") from e

    for col in ("uniprot", "reactant_smiles"):
        if col not in df.columns:
            raise DatasetError(
                f"в {csv_path.name} нет колонки '{col}'. Есть: {list(df.columns)}"
            )
    uniprot = str(df.iloc[0]["uniprot"]).strip()
    smiles = str(df.iloc[0]["reactant_smiles"]).strip()
    if not uniprot or not smiles:
        raise DatasetError("пустые uniprot/smiles в первой строке")
    return uniprot, smiles


def ligand_descriptors(smiles: str) -> tuple[Chem.Mol, float, float]:
    """SMILES → (mol, TPSA, LogP). RDKit-исключения оборачиваем в LigandError."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise LigandError(f"RDKit не распарсил SMILES: {smiles!r}")
    try:
        tpsa = Descriptors.TPSA(mol)
        logp = Descriptors.MolLogP(mol)
    except Exception as e:
        raise LigandError(f"RDKit не смог посчитать дескрипторы: {e}") from e
    logger.info("[ligand] TPSA=%.3f LogP=%.3f", tpsa, logp)
    return mol, tpsa, logp


def build_ligand_3d(mol: Chem.Mol) -> torch.Tensor:
    """SMILES → [N_atoms, 3] float32 через ETKDG."""
    mol_hs = Chem.AddHs(mol)
    try:
        rc = AllChem.EmbedMolecule(mol_hs, AllChem.ETKDG())
    except Exception as e:
        raise LigandError(f"EmbedMolecule упал: {e}") from e
    if rc != 0:
        logger.warning(
            "[ligand] ETKDG не сошёлся (rc=%d), пробую randomCoords=True", rc
        )
        rc = AllChem.EmbedMolecule(mol_hs, AllChem.ETKDG(), randomSeed=42)
        if rc != 0:
            raise LigandError("не удалось построить 3D-конформер лиганда")

    conformer = mol_hs.GetConformer()
    coords = [list(conformer.GetAtomPosition(i)) for i in range(mol_hs.GetNumAtoms())]
    tensor = torch.tensor(coords, dtype=torch.float32)
    logger.info("[ligand] 3D-граф субстрата: shape=%s", tuple(tensor.shape))
    return tensor


# =========================================================
# 6. Основной пайплайн
# =========================================================
def run_pipeline() -> int:
    """Возвращает exit-code: 0 = success, 1 = fallback, 2 = fatal."""
    t0 = time.time()
    log_environment()

    # --- 6.1 Датасет ---
    logger.info("Шаг 1/4: чтение CatPred-DB")
    try:
        uniprot, smiles = load_first_row(DATASET_PATH)
    except DatasetError as e:
        logger.error("❌ %s", e)
        return 2
    logger.info("[main] uniprot=%s smiles=%s", uniprot, smiles)

    # --- 6.2 Лиганд ---
    logger.info("Шаг 2/4: расчёт дескрипторов лиганда (RDKit)")
    try:
        mol, tpsa, logp = ligand_descriptors(smiles)
    except LigandError as e:
        logger.error("❌ %s", e)
        return 2

    # --- 6.3 Структура фермента ---
    logger.info("Шаг 3/4: загрузка 3D-структуры фермента")
    pdb_path = PDB_DIR / f"{uniprot}.pdb"
    model = PhysicsAdd3DGraphNet()

    coords: torch.Tensor | None = None
    source_label: str
    exit_code = 0

    try:
        logger.info("[main] ищу файл: %s", pdb_path.resolve())
        coords = model.extract_protein_coords(pdb_path)
        source_label = "РЕАЛЬНЫЕ 3D КООРДИНАТЫ ФЕРМЕНТА (AlphaFold)"
        logger.info("✅ использую структуру фермента")
    except ProteinFormatError as e:
        logger.error("❌ структура фермента невалидна: %s", e)
        logger.warning("⚠️  перехожу на бейзлайн: 3D-граф субстрата")
        exit_code = 1
    except ProteinParseError as e:
        logger.error("❌ парсинг PDB не удался: %s", e)
        logger.warning("⚠️  перехожу на бейзлайн: 3D-граф субстрата")
        exit_code = 1
    except FileNotFoundError as e:
        logger.error("❌ файл не найден: %s", e)
        logger.warning("⚠️  перехожу на бейзлайн: 3D-граф субстрата")
        exit_code = 1

    if coords is None or coords.numel() == 0:
        try:
            coords = build_ligand_3d(mol)
            source_label = "3D-ГРАФ СУБСТРАТА (fallback, без фермента)"
        except LigandError as e:
            logger.error("❌ fallback тоже упал: %s", e)
            return 2

    # --- 6.4 Forward pass ---
    logger.info("Шаг 4/4: forward-pass нейросети")
    try:
        pose_prob, center, barrier = model(coords, tpsa, logp)
    except Exception as e:
        logger.error("❌ forward упал: %s", e)
        logger.debug("traceback:\n%s", traceback.format_exc())
        return 2

    # --- 6.5 Отчёт ---
    elapsed = time.time() - t0
    print()
    print("=" * 62)
    print("📊 ФИЗИКО-ХИМИЧЕСКИЙ ОТЧЁТ СИСТЕМЫ")
    print("=" * 62)
    print(f"  Источник 3D-геометрии : {source_label}")
    print(f"  Субстрат (SMILES)     : {smiles}")
    print(f"  Фермент (UniProt)     : {uniprot}")
    print(
        f"  Файл структуры        : {pdb_path if exit_code == 0 else '— (не использован)'}"
    )
    print(f"  Количество 3D-узлов   : {coords.shape[0]}")
    print(f"  TPSA / LogP           : {tpsa:.2f} / {logp:.2f}")
    print(
        f"  Центр масс (X, Y, Z)  : [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}]"
    )
    print(f"  Барьер сольват. шубы  : {barrier:.2f} кДж/моль")
    print(f"  Вероятность посадки   : {pose_prob.item() * 100:.2f} %")
    print(f"  Время выполнения      : {elapsed:.2f} с")
    print("=" * 62)
    if exit_code == 1:
        print(
            "⚠️  ВНИМАНИЕ: работали на бейзлайне, а не на реальной структуре фермента."
        )
        print("    Проверь outputs/logs/main.log — там причина (HTML? DNS? VPN?).")
        print("=" * 62)
    print()
    logger.info("✓ пайплайн завершён за %.2f с, exit_code=%d", elapsed, exit_code)
    return exit_code


# =========================================================
# 7. Точка входа с глобальным обработчиком исключений
# =========================================================
def main() -> int:
    logger.info("=== СТАРТ Physics-Informed 3D Graph ML ===")
    try:
        return run_pipeline()
    except KeyboardInterrupt:
        logger.warning("Прервано пользователем (Ctrl+C)")
        return 130
    except Exception as e:
        logger.critical("НЕОЖИДАННОЕ ИСКЛЮЧЕНИЕ: %s", e)
        logger.critical("Traceback:\n%s", traceback.format_exc())
        logger.critical("Полный лог: %s", LOG_FILE.resolve())
        return 2


if __name__ == "__main__":
    sys.exit(main())
