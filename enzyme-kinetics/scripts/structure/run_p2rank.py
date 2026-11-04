from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import subprocess
from multiprocessing import Pool
import glob

BASE_DIR = str(PROJECT_DIR)
PDB_DIR = os.path.join(BASE_DIR, "data", "alphafold_pdb")  # Убедись, что путь верный!
OUT_DIR = os.path.join(BASE_DIR, "data", "processed", "p2rank_out")
PRANK_PATH = os.path.join(
    BASE_DIR, "vendor", "p2rank", "prank"
)  # Путь к исполняемому файлу


def run_prank(pdb_file):
    pdb_name = os.path.basename(pdb_file)
    out_file = os.path.join(OUT_DIR, pdb_name.replace(".pdb", ""))

    # Если уже обработано, пропускаем
    if os.path.exists(os.path.join(out_file, pdb_name + "_predictions.csv")):
        return f"Skip {pdb_name}"

    cmd = [
        PRANK_PATH,
        "predict",
        "-f",
        pdb_file,
        "-o",
        out_file,
        "-threads",
        "1",  # 1 поток на 1 файл, так как параллелим процессами
    ]
    try:
        subprocess.run(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
        )
        return f"Done {pdb_name}"
    except Exception:
        return f"Error {pdb_name}"


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    pdb_files = glob.glob(os.path.join(PDB_DIR, "*.pdb"))
    print(f"Найдено PDB файлов: {len(pdb_files)}")

    # Используем все ядра, кроме одного
    num_cores = min(4, os.cpu_count() or 1)
    print(f"Запуск P2Rank на {num_cores} ядрах...")

    with Pool(num_cores) as pool:
        for result in pool.imap_unordered(run_prank, pdb_files):
            print(result)
