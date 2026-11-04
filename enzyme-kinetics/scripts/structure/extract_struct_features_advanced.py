from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import pandas as pd
import re
from multiprocessing import Pool
import glob

BASE_DIR = str(PROJECT_DIR)
P2RANK_DIR = os.path.join(BASE_DIR, "data", "processed", "p2rank_out")
PDB_DIR = os.path.join(BASE_DIR, "data", "alphafold_pdb")
OUT_CSV = os.path.join(BASE_DIR, "data", "processed", "struct_features_advanced.csv")

CATALYTIC_AAS = ["HIS", "ASP", "GLU", "SER", "CYS", "TYR", "LYS", "ARG"]
POSITIVE_AAS = ["LYS", "ARG", "HIS"]
NEGATIVE_AAS = ["ASP", "GLU"]
HYDROPHOBIC_AAS = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]


def get_aas_from_pocket(pdb_path, target_resids):
    found_aas = []
    if not os.path.exists(pdb_path) or not target_resids:
        return found_aas
    seen = set()
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    resname = line[17:20].strip()
                    chain = line[21]
                    resid = line[22:26].strip()
                    res_key = f"{chain}_{resid}"
                    if res_key in target_resids and resname not in ("HOH", "WAT"):
                        if res_key not in seen:
                            found_aas.append(resname)
                            seen.add(res_key)
    except Exception:
        pass
    return found_aas


def process_single_pdb(csv_path):
    try:
        df = pd.read_csv(csv_path)
        df.columns = (
            df.columns.str.strip()
        )  # <-- фикс №1: убираем пробелы в именах колонок
        if df.empty:
            return None
        df = df.sort_values("rank")
        top_pocket = df.iloc[0]
        if pd.isna(top_pocket["residue_ids"]):
            return None
        score = float(top_pocket["score"])
        target_resids = set(str(top_pocket["residue_ids"]).split())
        pocket_size = len(target_resids)
        if pocket_size == 0:
            return None

        # <-- фикс №2: base с .pdb для поиска файла, pdb_name без .pdb для merge
        base = os.path.basename(csv_path).replace("_predictions.csv", "")
        pdb_name = base[:-4] if base.endswith(".pdb") else base
        pdb_path = os.path.join(PDB_DIR, base)
        if not os.path.exists(pdb_path):
            match = re.match(r"^AF-([A-Z0-9]+)-F1-model_v\d+$", pdb_name)
            if match:
                pdb_path = os.path.join(PDB_DIR, match.group(1) + ".pdb")
        if not os.path.exists(pdb_path):
            found = glob.glob(os.path.join(PDB_DIR, "**", base), recursive=True)
            pdb_path = found[0] if found else ""

        pocket_aas = get_aas_from_pocket(pdb_path, target_resids)
        total_aas = len(pocket_aas) if pocket_aas else 1

        features = {
            "pdb_name": pdb_name,
            "pocket_score": score,
            "pocket_size": pocket_size,
        }
        for aa in CATALYTIC_AAS:
            features[f"freq_{aa}"] = (
                pocket_aas.count(aa) / total_aas if pocket_aas else float("nan")
            )
        features["pocket_charge"] = (
            sum(pocket_aas.count(aa) for aa in POSITIVE_AAS)
            - sum(pocket_aas.count(aa) for aa in NEGATIVE_AAS)
        ) / total_aas
        features["pocket_hydrophobicity"] = (
            sum(pocket_aas.count(aa) for aa in HYDROPHOBIC_AAS) / total_aas
        )
        if not pocket_aas:
            features["pocket_charge"] = float("nan")
            features["pocket_hydrophobicity"] = float("nan")
        return features
    except Exception:
        return None


if __name__ == "__main__":
    csv_files = glob.glob(
        os.path.join(P2RANK_DIR, "**", "*_predictions.csv"), recursive=True
    )
    print(f"Найдено результатов P2Rank: {len(csv_files)}")

    num_cores = min(4, os.cpu_count() or 1)
    print(f"Запускаем парсинг PDB на {num_cores} ядрах...")

    results = []
    with Pool(num_cores) as pool:
        for i, res in enumerate(pool.imap_unordered(process_single_pdb, csv_files)):
            if res:
                results.append(res)
            if (i + 1) % 1000 == 0:
                print(f"Обработано {i + 1} / {len(csv_files)}")

    df_features = pd.DataFrame(results)
    if df_features.empty:
        print("⚠️ Ни одной строки не получилось — проверь пути P2RANK_DIR / PDB_DIR.")
    else:
        expected_cols = (
            ["pdb_name", "pocket_score", "pocket_size"]
            + [f"freq_{aa}" for aa in CATALYTIC_AAS]
            + ["pocket_charge", "pocket_hydrophobicity"]
        )
        for col in expected_cols:
            if col not in df_features.columns:
                df_features[col] = 0.0
        df_features = df_features[expected_cols]
        df_features.to_csv(OUT_CSV, index=False)
        print(f"✅ Готово! Сохранено {len(df_features)} строк в {OUT_CSV}")
        print(
            f"Размерность признаков: {df_features.shape[1] - 1} (без учёта имени PDB)"
        )
