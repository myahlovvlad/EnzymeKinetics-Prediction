from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed


def main():
    BASE = str(PROJECT_DIR)
    CIF = os.path.join(BASE, "data", "alphafold_cif")
    os.makedirs(CIF, exist_ok=True)

    with open(os.path.join(BASE, "data", "manifests", "missing_uids.txt")) as f:
        to_download = [line.strip() for line in f if line.strip()]

    print(f"Нужно скачать CIF: {len(to_download)}")

    HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    stats = {"ok": 0, "missing": 0, "fail": 0}

    def download(uid):
        url = f"https://storage.googleapis.com/public-datasets-deepmind-alphafold-v4/AF-{uid}-F1-model_v4.cif"
        path = os.path.join(CIF, f"{uid}.cif")
        try:
            r = requests.get(url, headers=HEADERS, stream=True, timeout=30)
            if r.status_code == 200:
                with open(path, "wb") as f:
                    for chunk in r.iter_content(8192):
                        f.write(chunk)
                return "ok"
            if r.status_code == 404:
                return "missing"
            return "fail"
        except Exception:
            return "fail"

    with ThreadPoolExecutor(max_workers=16) as ex:
        futures = {ex.submit(download, u): u for u in to_download}
        for i, fut in enumerate(as_completed(futures), 1):
            stats[fut.result()] += 1
            if i % 100 == 0 or i == len(to_download):
                print(
                    f"[{i}/{len(to_download)}] ok={stats['ok']} missing={stats['missing']} fail={stats['fail']}"
                )

    print(f"\n✅ ok={stats['ok']} missing={stats['missing']} fail={stats['fail']}")


if __name__ == "__main__":
    main()
