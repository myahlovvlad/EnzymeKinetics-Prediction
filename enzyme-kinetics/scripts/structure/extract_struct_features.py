from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import pandas as pd
import glob

BASE_DIR = str(PROJECT_DIR)
P2RANK_DIR = os.path.join(BASE_DIR, "data", "processed", "p2rank_out")
OUT_CSV = os.path.join(BASE_DIR, "data", "processed", "struct_features.csv")

# Список каталитических аминокислот (можно расширить)
CATALYTIC_AAS = ["HIS", "ASP", "GLU", "SER", "CYS", "TYR", "LYS", "ARG"]


def extract_features_from_p2rank(csv_path):
    try:
        # P2Rank выдает CSV с колонками: name, rank, score, probability, ..., residue_ids
        df = pd.read_csv(csv_path)
        if df.empty:
            return None

        # Берем топ-1 карман
        top_pocket = df.iloc[0]
        score = top_pocket["score"]

        # residue_ids выглядят как "A_45 A_46 A_78 ..."
        res_ids = str(top_pocket["residue_ids"]).split()

        # Извлекаем трехбуквенные коды аминокислот
        # Пример: "A_45" -> нам нужно понять, какая это аминокислота.
        # Но P2Rank выдает только номера, а не названия!
        # Поэтому мы вернем просто количество остатков в кармане и score.
        # (Для полного разбора AA нам нужно будет спарсить сам PDB, но начнем с простого)

        pocket_size = len(res_ids)

        return {
            "pdb_name": os.path.basename(csv_path).replace("_predictions.csv", ""),
            "pocket_score": score,
            "pocket_size": pocket_size,
        }
    except Exception:
        return None


if __name__ == "__main__":
    csv_files = glob.glob(
        os.path.join(P2RANK_DIR, "**", "*_predictions.csv"), recursive=True
    )
    print(f"Найдено результатов P2Rank: {len(csv_files)}")

    features = []
    for f in csv_files:
        res = extract_features_from_p2rank(f)
        if res:
            features.append(res)

    df_features = pd.DataFrame(features)
    df_features.to_csv(OUT_CSV, index=False)
    print(f"Сохранено в {OUT_CSV}. Размерность: {df_features.shape}")
