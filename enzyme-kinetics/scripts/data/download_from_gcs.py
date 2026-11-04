from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed


def main():
    BASE_DIR = str(PROJECT_DIR)
    DATA_DIR = os.path.join(BASE_DIR, "data", "Baseline", "kcat")
    PDB_DIR = os.path.join(BASE_DIR, "data", "alphafold_pdb")
    os.makedirs(PDB_DIR, exist_ok=True)

    print("Собираем список UniProt ID...")
    train = pd.read_csv(os.path.join(DATA_DIR, "kcat-random_trainval.csv"))
    test = pd.read_csv(os.path.join(DATA_DIR, "kcat-random_test.csv"))
    all_uniprots = pd.concat([train["uniprot"], test["uniprot"]]).dropna().unique()

    to_download = []
    for uid in all_uniprots:
        uid = str(uid).strip()
        if not os.path.exists(os.path.join(PDB_DIR, f"{uid}.pdb")):
            to_download.append(uid)

    print(f"Нужно скачать: {len(to_download)} файлов.")

    HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    stats = {"ok": 0, "missing": 0, "fail": 0}

    def download_single(uid):
        url = f"https://storage.googleapis.com/alphafold-databases/v4/AF-{uid}-F1-model_v4.pdb"
        filepath = os.path.join(PDB_DIR, f"{uid}.pdb")
        try:
            r = requests.get(url, headers=HEADERS, stream=True, timeout=30)
            if r.status_code == 200:
                with open(filepath, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                return "ok"
            elif r.status_code == 404:
                return "missing"
            else:
                return "fail"
        except Exception:
            return "fail"

    # 16 потоков — Google это спокойно тянет
    with ThreadPoolExecutor(max_workers=16) as ex:
        futures = {ex.submit(download_single, uid): uid for uid in to_download}
        for i, fut in enumerate(as_completed(futures), 1):
            stats[fut.result()] += 1
            if i % 100 == 0 or i == len(to_download):
                print(
                    f"[{i}/{len(to_download)}] ok={stats['ok']} missing={stats['missing']} fail={stats['fail']}"
                )

    print(
        f"\n✅ Готово. Скачано: {stats['ok']}, отсутствует в базе: {stats['missing']}, ошибок: {stats['fail']}"
    )


if __name__ == "__main__":
    main()
