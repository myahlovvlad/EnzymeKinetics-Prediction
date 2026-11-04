#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Аудит датасета enzyme-kinetics (kcat / km / ki).
Запуск в WSL (через uv):
    cd /home/vlad_myahlov/projects/enzyme-kinetics
    uv run python audit_enzyme_data.py            # полный аудит
    uv run python audit_enzyme_data.py --fast     # без чтения PT/PKL больше 500 MB

Ничего не модифицирует в data/ — только читает.
Отчёт пишется в ./audit_out/
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
import sys
import traceback
from datetime import datetime
from pathlib import Path

from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR
from typing import Any, Optional

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Конфиг
# --------------------------------------------------------------------------- #

DATA_ROOT = DATA_DIR
OUT_DIR = OUTPUT_DIR / "audit"

TARGETS = ["kcat", "km", "ki"]

SEQUENCE_COL = "sequence"
PDB_COL = "pdbpath"
UNIPROT_COL = "uniprot"
SMILES_COL = "reaction_smiles"  # только у kcat
SUBSTRATE_COL = "substrate_smiles"  # только у km / ki

TARGET_COLS = {
    "kcat": ["value", "log10_value", "log10kcat_max"],
    "km": ["value", "log10_value", "log10km_mean"],
    "ki": ["value", "log10_value", "log10ki_mean"],
}
# Это выходы модели/бейзлайна — как признаки ИСПОЛЬЗОВАТЬ НЕЛЬЗЯ
PREDICTED_COLS = {
    "kcat": ["predicted_log10kcat_max"],
    "km": ["predicted_log10km_mean"],
    "ki": ["predicted_log10ki_mean"],
}
CLUSTER_SUFFIXES = ["20cluster", "40cluster", "60cluster", "80cluster", "99cluster"]
GROUP_COLS = ["group", "top_3_seqs", "top_3_seqids"]


def csv_files_for(target: str, data_root: Optional[Path] = None) -> dict:
    base = (data_root or DATA_ROOT) / "Baseline" / target
    return {
        "trainval": base / f"{target}-random_trainval.csv",
        "test": base / f"{target}-random_test.csv",
        "test_mean3sim": base / f"{target}-random_test_mean_3sim.csv",
        **{
            f"test_{sfx}_mean3sim": base
            / f"{target}-random_test_sequence_{sfx}_mean_3sim.csv"
            for sfx in CLUSTER_SUFFIXES
        },
        **{
            f"seq_test_{sfx}": base / f"{target}-seq_test_sequence_{sfx}.csv"
            for sfx in CLUSTER_SUFFIXES
        },
        "identities_pkl": base / f"{target}_test_train_identities_updated.pkl",
    }


def build_extra_files(data_root: Path) -> dict:
    return {
        "catpred_pt_top": data_root / "CatPred-DB" / "catpred_progres_embeds_dict.pt",
        "catpred_pt_data": data_root
        / "CatPred-DB"
        / "data"
        / "catpred_progres_embeds_dict.pt",
        "test_seqset_json": data_root / "CatPred-DB" / "data" / "test_seqset.json",
        "single_pdb": data_root / "alphafold_pdb" / "Q05469.pdb",
    }


# Значение по умолчанию; в main() пересобирается при --data-root
EXTRA_FILES = build_extra_files(DATA_ROOT)


# --------------------------------------------------------------------------- #
# Логирование
# --------------------------------------------------------------------------- #


def setup_logging() -> logging.Logger:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fmt = "%(asctime)s | %(levelname)-7s | %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(OUT_DIR / "audit.log", mode="w", encoding="utf-8"),
        ],
    )
    return logging.getLogger("audit")


def log_environment(logger: logging.Logger) -> None:
    logger.info(f"python = {sys.version.split()[0]}  ({sys.executable})")
    logger.info(f"uv-run = {'UV_RUN_RECURSION_DEPTH' in os.environ}")
    logger.info(f"pandas = {pd.__version__}")
    logger.info(f"numpy  = {np.__version__}")
    try:
        import torch

        logger.info(f"torch  = {torch.__version__}  cuda={torch.cuda.is_available()}")
    except ImportError:
        logger.info("torch  = not installed (only CSV/JSON/PKL audit)")


def section(logger: logging.Logger, title: str) -> None:
    logger.info("")
    logger.info("=" * 88)
    logger.info(title)
    logger.info("=" * 88)


def fmt_size(n: int) -> str:
    x = float(n)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if x < 1024:
            return f"{x:.2f} {unit}"
        x /= 1024
    return f"{x:.2f} PB"


# --------------------------------------------------------------------------- #
# Утилиты
# --------------------------------------------------------------------------- #


def safe_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def file_meta(path: Path) -> dict:
    if not path.exists():
        return {"path": str(path), "exists": False}
    st = path.stat()
    return {
        "path": str(path),
        "exists": True,
        "size_bytes": st.st_size,
        "size_human": fmt_size(st.st_size),
        "mtime": datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
    }


def resolve_pdb_path(pdb_str: Any, data_root: Path) -> Optional[Path]:
    """
    Пробуем разрешить pdbpath из CSV в реальный файл.
    CSV содержит варианты:
        ./data/raw/alphafold//AF-XXXX-F1-model_v4.pdb
        ./data/raw/esmfold//ESMFold-XXXX-v1.pdb
    А в data/ ещё есть папка alphafold_pdb/.
    """
    if pdb_str is None:
        return None
    if isinstance(pdb_str, float) and np.isnan(pdb_str):
        return None
    s = str(pdb_str).strip()
    if not s:
        return None
    s = s.replace("\\", "/")

    candidates: list[Path] = []
    p = Path(s)
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.append(data_root / s.lstrip("./"))
        if s.startswith("./data/"):
            candidates.append(data_root / s[len("./data/") :])
        if s.startswith("data/"):
            candidates.append(data_root / s[len("data/") :])
        candidates.append(data_root.parent / s.lstrip("./"))

    basename = Path(s).name
    for sub in ["raw/alphafold", "raw/esmfold", "alphafold_pdb", "raw", ""]:
        candidates.append(data_root / sub / basename if sub else data_root / basename)

    for c in candidates:
        if c.exists():
            return c
    return None


# --------------------------------------------------------------------------- #
# CSV audit
# --------------------------------------------------------------------------- #


def audit_csv(logger: logging.Logger, path: Path, target: str, kind: str) -> dict:
    res: dict[str, Any] = {"kind": kind, "target": target, **file_meta(path)}
    if not path.exists():
        logger.warning(f"  [MISSING] {path}")
        res["status"] = "missing"
        return res

    logger.info(f"  {kind:24s} {path.name}  ({res['size_human']})")
    try:
        df = pd.read_csv(path, low_memory=False)
    except Exception as e:
        logger.error(f"    ! read_csv failed: {e}")
        res["status"] = "read_error"
        res["error"] = str(e)
        return res

    res["status"] = "ok"
    res["n_rows"] = int(len(df))
    res["n_cols"] = int(df.shape[1])
    res["columns"] = list(df.columns)
    logger.info(f"    shape: {df.shape}")

    required = [SEQUENCE_COL] + TARGET_COLS[target]
    missing_required = [c for c in required if c not in df.columns]
    if missing_required:
        logger.warning(f"    ! отсутствуют обязательные столбцы: {missing_required}")
    res["missing_required"] = missing_required

    nulls = {}
    for c in required + [PDB_COL]:
        if c in df.columns:
            nulls[c] = int(df[c].isna().sum())
    res["nulls"] = nulls
    logger.info(f"    nulls: {nulls}")

    if SEQUENCE_COL in df.columns:
        seqs = df[SEQUENCE_COL].astype(str)
        empty = int((seqs.str.strip() == "").sum())
        dups = int(df[SEQUENCE_COL].duplicated().sum())
        logger.info(
            f"    seq: total={len(df)}, unique={df[SEQUENCE_COL].nunique()}, "
            f"dups={dups}, empty={empty}"
        )
        res["dup_sequences"] = dups
        res["empty_sequences"] = empty
        res["unique_sequences"] = int(df[SEQUENCE_COL].nunique())

        aa_regex_ok = seqs.str.fullmatch(r"[A-Za-z*Xx\-\.]*")
        bad = int((~aa_regex_ok).sum())
        if bad:
            logger.warning(f"    ! подозрительные символы в sequence у {bad} строк")
        res["bad_seq_chars"] = bad

    res["target_stats"] = {}
    for c in TARGET_COLS[target]:
        if c in df.columns:
            s = safe_numeric(df[c])
            n_ok = int(s.notna().sum())
            if n_ok > 0:
                stats = {
                    "n_nonnull": n_ok,
                    "mean": float(s.mean()),
                    "std": float(s.std()) if n_ok > 1 else 0.0,
                    "min": float(s.min()),
                    "max": float(s.max()),
                }
                logger.info(
                    f"    {c:16s} n={n_ok}, mean={stats['mean']:.3f}, "
                    f"std={stats['std']:.3f}, min={stats['min']}, max={stats['max']}"
                )
            else:
                stats = {"n_nonnull": 0}
                logger.info(f"    {c:16s} all-NaN")
            res["target_stats"][c] = stats

    present_pred = [c for c in PREDICTED_COLS[target] if c in df.columns]
    if present_pred:
        logger.info(f"    predicted cols (НЕ признаки): {present_pred}")
    res["predicted_cols_present"] = present_pred

    if "group" in df.columns:
        n_groups = int(df["group"].nunique())
        logger.info(f"    groups: {n_groups}")
        res["n_groups"] = n_groups
    if "top_3_seqs" in df.columns and df["top_3_seqs"].notna().any():
        sample = str(df["top_3_seqs"].dropna().iloc[0])
        res["top_3_sample_head"] = sample[:120]

    if PDB_COL in df.columns:
        sample = df[PDB_COL].dropna().astype(str).head(3).tolist()
        res["pdbpath_sample"] = sample
        resolved = 0
        unresolved_examples: list[str] = []
        for v in df[PDB_COL].dropna().astype(str):
            if resolve_pdb_path(v, DATA_ROOT) is not None:
                resolved += 1
            elif len(unresolved_examples) < 3:
                unresolved_examples.append(v)
        n_pdb = int(df[PDB_COL].notna().sum())
        logger.info(f"    pdbpath: {resolved}/{n_pdb} разрешено в реальный файл")
        if unresolved_examples:
            logger.warning(f"    ! примеры нерешённых pdbpath: {unresolved_examples}")
        res["pdb_resolved"] = resolved
        res["pdb_nonnull"] = n_pdb
        res["pdb_unresolved_examples"] = unresolved_examples

    return res


# --------------------------------------------------------------------------- #
# JSON / PKL / PT
# --------------------------------------------------------------------------- #


def audit_json(logger: logging.Logger, path: Path) -> dict:
    res = {"kind": "json", **file_meta(path)}
    if not path.exists():
        logger.warning(f"  [MISSING] {path}")
        res["status"] = "missing"
        return res
    logger.info(f"  json  {path.name}  ({res['size_human']})")
    try:
        with open(path, "r", encoding="utf-8") as f:
            obj = json.load(f)
        res["status"] = "ok"
        if isinstance(obj, dict):
            res["top_keys"] = list(obj.keys())
            res["per_key"] = {}
            for k, v in obj.items():
                if isinstance(v, list):
                    seqs = [str(x) for x in v]
                    logger.info(
                        f"    {k}: list, n={len(seqs)}, unique={len(set(seqs))}"
                    )
                    res["per_key"][k] = {
                        "type": "list",
                        "n": len(seqs),
                        "unique": len(set(seqs)),
                    }
                else:
                    logger.info(f"    {k}: {type(v).__name__}")
                    res["per_key"][k] = {"type": type(v).__name__}
        elif isinstance(obj, list):
            res["n"] = len(obj)
            logger.info(f"    list, n={len(obj)}")
        else:
            res["type"] = type(obj).__name__
            logger.info(f"    type={type(obj).__name__}")
    except Exception as e:
        logger.error(f"    ! json load failed: {e}")
        res["status"] = "error"
        res["error"] = str(e)
    return res


def audit_pickle(logger: logging.Logger, path: Path, fast: bool) -> dict:
    res = {"kind": "pickle", **file_meta(path)}
    if not path.exists():
        logger.warning(f"  [MISSING] {path}")
        res["status"] = "missing"
        return res
    size = res["size_bytes"]
    logger.info(f"  pkl   {path.name}  ({res['size_human']})")
    if fast and size > 500 * 1024 * 1024:
        logger.warning("    ! --fast и размер > 500 MB — пропускаю загрузку")
        res["status"] = "skipped_size_guard"
        return res
    try:
        with open(path, "rb") as f:
            obj = pickle.load(f)
        res["status"] = "ok"
        res["python_type"] = type(obj).__name__
        logger.info(f"    type={type(obj).__name__}")
        if isinstance(obj, dict):
            keys = list(obj.keys())[:10]
            logger.info(f"    dict, size={len(obj)}, sample keys={keys}")
            res["dict_size"] = len(obj)
            res["sample_keys"] = [str(k) for k in keys]
            if obj:
                first_val = next(iter(obj.values()))
                logger.info(f"    first value type: {type(first_val).__name__}")
                res["first_value_type"] = type(first_val).__name__
        elif isinstance(obj, pd.DataFrame):
            logger.info(
                f"    DataFrame shape={obj.shape}, cols={list(obj.columns)[:10]}"
            )
            res["dataframe_shape"] = list(obj.shape)
            res["dataframe_columns"] = list(obj.columns)
        elif isinstance(obj, np.ndarray):
            logger.info(f"    ndarray shape={obj.shape}, dtype={obj.dtype}")
            res["ndarray_shape"] = list(obj.shape)
            res["ndarray_dtype"] = str(obj.dtype)
        elif isinstance(obj, (list, tuple)):
            logger.info(f"    {type(obj).__name__} len={len(obj)}")
            res["seq_len"] = len(obj)
    except Exception as e:
        logger.error(f"    ! pickle load failed: {e}")
        res["status"] = "error"
        res["error"] = str(e)
    return res


def audit_torch(logger: logging.Logger, path: Path, fast: bool) -> dict:
    res = {"kind": "torch_pt", **file_meta(path)}
    if not path.exists():
        logger.warning(f"  [MISSING] {path}")
        res["status"] = "missing"
        return res
    size = res["size_bytes"]
    logger.info(f"  pt    {path.name}  ({res['size_human']})")
    try:
        import torch
    except ImportError:
        logger.warning(
            "    ! torch не установлен — только размер. "
            "Установите: uv add torch --index-url https://download.pytorch.org/whl/cpu"
        )
        res["status"] = "torch_not_installed"
        return res

    if fast and size > 500 * 1024 * 1024:
        logger.warning("    ! --fast и размер > 500 MB — пропускаю загрузку")
        res["status"] = "skipped_size_guard"
        return res

    try:
        try:
            obj = torch.load(path, map_location="cpu", weights_only=False)
        except TypeError:
            obj = torch.load(path, map_location="cpu")
        res["status"] = "ok"
        res["python_type"] = type(obj).__name__
        logger.info(f"    type={type(obj).__name__}")
        if isinstance(obj, dict):
            n = len(obj)
            sample_keys = list(obj.keys())[:5]
            logger.info(f"    dict, size={n}, sample keys={sample_keys}")
            res["dict_size"] = n
            res["sample_keys"] = [str(k) for k in sample_keys]
            if n:
                v = next(iter(obj.values()))
                if hasattr(v, "shape"):
                    logger.info(
                        f"    first value shape={tuple(v.shape)}, dtype={v.dtype}"
                    )
                    res["first_value_shape"] = list(v.shape)
                    res["first_value_dtype"] = str(v.dtype)
                elif isinstance(v, (list, tuple)):
                    logger.info(f"    first value {type(v).__name__}, len={len(v)}")
                    res["first_value_len"] = len(v)
                else:
                    logger.info(f"    first value type={type(v).__name__}")
        elif hasattr(obj, "shape"):
            logger.info(f"    tensor shape={tuple(obj.shape)}, dtype={obj.dtype}")
            res["tensor_shape"] = list(obj.shape)
            res["tensor_dtype"] = str(obj.dtype)
    except Exception as e:
        logger.error(f"    ! torch.load failed: {e}")
        res["status"] = "error"
        res["error"] = str(e)
    return res


# --------------------------------------------------------------------------- #
# Overlap analysis
# --------------------------------------------------------------------------- #


def overlap_sets(a: pd.Series, b: pd.Series) -> dict:
    sa = set(a.dropna().astype(str))
    sb = set(b.dropna().astype(str))
    inter = sa & sb
    return {
        "n_a": len(sa),
        "n_b": len(sb),
        "n_overlap": len(inter),
        "jaccard": (len(inter) / len(sa | sb)) if (sa | sb) else 0.0,
        "examples": list(inter)[:3],
    }


def run_overlap_analysis(
    logger: logging.Logger, dfs: dict[str, pd.DataFrame], target: str
) -> dict:
    res: dict[str, Any] = {}
    logger.info(f"--- overlaps for {target} ---")
    train = dfs.get("trainval")
    if train is None:
        logger.warning("  trainval отсутствует, overlap пропущен")
        return res

    kinds = (
        ["test", "test_mean3sim"]
        + [f"test_{s}_mean3sim" for s in CLUSTER_SUFFIXES]
        + [f"seq_test_{s}" for s in CLUSTER_SUFFIXES]
    )

    for kind in kinds:
        test = dfs.get(kind)
        if test is None:
            continue
        o = overlap_sets(train[SEQUENCE_COL], test[SEQUENCE_COL])
        logger.info(
            f"  trainval ∩ {kind:28s} : "
            f"|train|={o['n_a']}, |{kind}|={o['n_b']}, "
            f"overlap={o['n_overlap']}, jaccard={o['jaccard']:.4f}"
        )
        res[kind] = o

        if UNIPROT_COL in train.columns and UNIPROT_COL in test.columns:
            ou = overlap_sets(train[UNIPROT_COL], test[UNIPROT_COL])
            logger.info(
                f"    uniprot overlap: {ou['n_overlap']} (jaccard={ou['jaccard']:.4f})"
            )
            res[kind]["uniprot_overlap"] = ou["n_overlap"]

    return res


# --------------------------------------------------------------------------- #
# Дерево директорий
# --------------------------------------------------------------------------- #


def print_tree(logger: logging.Logger) -> None:
    section(logger, "Directory tree (top-2 уровня + счётчики)")
    if not DATA_ROOT.exists():
        logger.error(f"DATA_ROOT не существует: {DATA_ROOT}")
        return
    for p in sorted(DATA_ROOT.iterdir()):
        if p.is_dir():
            files = [f for f in p.rglob("*") if f.is_file()]
            n_pdb = sum(1 for f in files if f.suffix.lower() == ".pdb")
            total = sum(f.stat().st_size for f in files)
            logger.info(
                f"  {p.name}/  files={len(files)}  "
                f"(pdb={n_pdb}, size={fmt_size(total)})"
            )
            if len(files) <= 12:
                for f in sorted(files):
                    logger.info(
                        f"    - {f.relative_to(DATA_ROOT)}  "
                        f"({fmt_size(f.stat().st_size)})"
                    )
        else:
            logger.info(f"  {p.name}  ({fmt_size(p.stat().st_size)})")


# --------------------------------------------------------------------------- #
# Итоговый pipeline
# --------------------------------------------------------------------------- #


def print_pipeline_recommendation(logger: logging.Logger) -> None:
    section(logger, "РЕКОМЕНДУЕМЫЙ ПАЙПЛАЙН ЗАГРУЗКИ (без утечек)")
    logger.info("""
1) Единая точка входа
   DATA_ROOT = /home/vlad_myahlov/projects/enzyme-kinetics/data
   Для каждого target ∈ {kcat, km, ki} — CSV из data/Baseline/<target>/.

2) trainval / test
   train = <target>-random_trainval.csv
   test  = <target>-random_test.csv
   Cluster-файлы НЕ подмешивать в train/val — это отдельные тест-срезы.

3) Целевая переменная
   kcat  -> log10kcat_max
   km    -> log10km_mean
   ki    -> log10ki_mean
   log10_value — лог отдельного измерения; для CatPred-бенчмарка цель —
   агрегированное max/mean.
   Строки с NaN в цели — drop ТОЛЬКО для соответствующего target.

4) Запрещённые признаки (data leakage)
   predicted_log10kcat_max / predicted_log10km_mean / predicted_log10ki_mean
   group, top_3_seqs, top_3_seqids — служебные.
   *_cluster столбцы — ID соседей/кластеров для split, не признаки.

5) Валидация без утечки
   random test — базовый бенчмарк (пересечения последовательностей ожидаемы).
   Честная оценка обобщения — на cluster-срезах:
      20cluster -> строгая (почти нет гомологов)
      99cluster -> мягкая
   Отчёт по всем: random + 20/40/60/80/99.

6) PDB и структуры
   resolve_pdb_path() -> если файл существует, можно кодировать структуру
   (или брать эмбеддинги из catpred_progres_embeds_dict.pt по ключу =
   последовательность).
   Если PDB нет — не «нули молча», а маска has_structure=0 и/или отдельная
   ветка модели.

7) Эмбеддинги из .pt
   catpred_progres_embeds_dict.pt: dict sequence -> tensor.
   torch.load(..., map_location="cpu", weights_only=False).
   Соответствие с CSV проверять по строке sequence 1:1.

8) Дедупликация
   Внутри train удалить дубликаты по sequence (или по (sequence, reaction/
   substrate_smiles), если задача per reaction).
   Между train и test при random split — оставить как есть, но обязательно
   залогировать величину пересечения.

9) splits
   fit  <- trainval
   val  <- 5–10% из trainval (стратифицированно по target-бинам)
   test <- <target>-random_test.csv + <target>-random_test_sequence_*cluster*.csv
   Никаких пересечений между fit и val по sequence.

10) Отчётность
    MAE и R² по log10-цели, отдельно по random и по каждому cluster-срезу.
    Для km/ki сообщать в log10(M) — не переводить обратно без оговорки.
""")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def main() -> int:
    global DATA_ROOT, EXTRA_FILES

    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="не грузить PT/PKL > 500 MB")
    ap.add_argument(
        "--data-root", type=Path, default=None, help=f"по умолчанию {DATA_ROOT}"
    )
    args = ap.parse_args()

    if args.data_root is not None:
        DATA_ROOT = args.data_root.expanduser().resolve()
    EXTRA_FILES = build_extra_files(DATA_ROOT)

    logger = setup_logging()
    logger.info(f"DATA_ROOT = {DATA_ROOT}")
    logger.info(f"OUT_DIR   = {OUT_DIR}")
    logger.info(f"started   = {datetime.now().isoformat(timespec='seconds')}")
    log_environment(logger)

    report: dict[str, Any] = {"data_root": str(DATA_ROOT), "targets": {}}

    print_tree(logger)

    for target in TARGETS:
        section(logger, f"TARGET = {target}")
        files = csv_files_for(target)
        target_report: dict[str, Any] = {"files": {}}
        dfs: dict[str, pd.DataFrame] = {}

        for kind, path in files.items():
            if kind == "identities_pkl":
                continue
            meta = audit_csv(logger, path, target, kind)
            target_report["files"][kind] = meta
            if meta.get("status") == "ok":
                try:
                    dfs[kind] = pd.read_csv(path, low_memory=False)
                except Exception:
                    pass

        section(logger, f"OVERLAP ANALYSIS — {target}")
        target_report["overlap"] = run_overlap_analysis(logger, dfs, target)

        section(logger, f"IDENTITIES PKL — {target}")
        target_report["identities_pkl"] = audit_pickle(
            logger, files["identities_pkl"], fast=args.fast
        )

        report["targets"][target] = target_report

    section(logger, "EXTRA FILES")
    report["extras"] = {}
    for name, path in EXTRA_FILES.items():
        if path.suffix == ".json":
            report["extras"][name] = audit_json(logger, path)
        elif path.suffix == ".pkl":
            report["extras"][name] = audit_pickle(logger, path, fast=args.fast)
        elif path.suffix == ".pt":
            report["extras"][name] = audit_torch(logger, path, fast=args.fast)
        elif path.suffix == ".pdb":
            m = file_meta(path)
            logger.info(f"  pdb   {path.name}  ({m.get('size_human', 'missing')})")
            report["extras"][name] = {"kind": "pdb", **m}
        else:
            m = file_meta(path)
            logger.info(f"  ?     {path.name}  ({m.get('size_human', 'missing')})")
            report["extras"][name] = {"kind": "unknown", **m}

    section(logger, "CROSS-CHECK: test_seqset.json vs CSV test")
    json_meta = report["extras"].get("test_seqset_json", {})
    if json_meta.get("status") == "ok" and json_meta.get("per_key"):
        try:
            with open(EXTRA_FILES["test_seqset_json"], "r", encoding="utf-8") as f:
                seqset = json.load(f)
        except Exception:
            seqset = {}
        cross: dict[str, Any] = {}
        for target in TARGETS:
            if target not in seqset or not isinstance(seqset[target], list):
                continue
            js = set(map(str, seqset[target]))
            test_path = csv_files_for(target)["test"]
            if not test_path.exists():
                continue
            tdf = pd.read_csv(test_path, low_memory=False)
            ts = set(tdf[SEQUENCE_COL].dropna().astype(str))
            inter = js & ts
            logger.info(
                f"  {target:5s}: json={len(js)}, csv_test={len(ts)}, "
                f"overlap={len(inter)}"
            )
            cross[target] = {
                "json": len(js),
                "csv_test": len(ts),
                "overlap": len(inter),
            }
        report["test_seqset_crosscheck"] = cross

    section(logger, "SUMMARY")
    for target, tr in report["targets"].items():
        train = tr["files"].get("trainval", {})
        test = tr["files"].get("test", {})
        logger.info(
            f"  {target:5s}: train_rows={train.get('n_rows', '?')}, "
            f"test_rows={test.get('n_rows', '?')}, "
            f"train_pdb_resolved={train.get('pdb_resolved', '?')}/"
            f"{train.get('pdb_nonnull', '?')}"
        )
        ov = tr["overlap"].get("test", {})
        if ov:
            logger.info(
                f"         trainval ∩ test (sequence): "
                f"{ov['n_overlap']} (jaccard={ov['jaccard']:.4f})"
            )

    report_path = OUT_DIR / "audit_summary.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2, default=str)
    logger.info(f"\nJSON-отчёт: {report_path}")
    logger.info(f"LOG:        {OUT_DIR / 'audit.log'}")

    try:
        import torch  # noqa: F401
    except ImportError:
        logger.warning(
            "torch не установлен — .pt-файлы не проверены. "
            "Установите: uv add torch --index-url https://download.pytorch.org/whl/cpu"
        )

    print_pipeline_recommendation(logger)
    logger.info("done.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
