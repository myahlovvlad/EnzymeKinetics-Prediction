from enzymekinetics_prediction.paths import PROJECT_DIR

# scripts/download_missing_ebi.py
import os
import time
import random
import requests


def main():
    BASE = str(PROJECT_DIR)
    PDB = os.path.join(BASE, "data", "alphafold_pdb")
    os.makedirs(PDB, exist_ok=True)

    with open(os.path.join(BASE, "data", "manifests", "missing_uids.txt")) as f:
        uids = [line.strip() for line in f if line.strip()]

    print(f"К скачиванию: {len(uids)} ID")

    HEADERS = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }
    stats = {"ok": 0, "exists": 0, "fail": 0, "404": 0}

    for i, uid in enumerate(uids, 1):
        path = os.path.join(PDB, f"{uid}.pdb")
        if os.path.exists(path):
            stats["exists"] += 1
            continue

        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v6.pdb"
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            if r.status_code == 200:
                with open(path, "wb") as f:
                    f.write(r.content)
                stats["ok"] += 1
            elif r.status_code == 404:
                stats["404"] += 1
            else:
                stats["fail"] += 1
        except Exception:
            stats["fail"] += 1

        if i % 50 == 0 or i == len(uids):
            print(
                f"[{i}/{len(uids)}] ok={stats['ok']} exists={stats['exists']} 404={stats['404']} fail={stats['fail']}"
            )

        # Случайная пауза 3–6 секунд
        time.sleep(random.uniform(3, 6))

    print(
        f"\n✅ Готово. ok={stats['ok']} exists={stats['exists']} 404={stats['404']} fail={stats['fail']}"
    )


if __name__ == "__main__":
    main()
