from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import re


def main():
    BASE = str(PROJECT_DIR)
    PDB = os.path.join(BASE, "data", "alphafold_pdb")

    # Паттерн: AF-{uid}-F1-model_v4.pdb  →  {uid}.pdb
    pat = re.compile(r"^AF-([A-Z0-9]+)-F1-model_v\d+\.pdb$", re.IGNORECASE)

    renamed, skipped = 0, 0
    for name in os.listdir(PDB):
        if not name.endswith(".pdb"):
            continue
        m = pat.match(name)
        if not m:
            skipped += 1
            continue
        uid = m.group(1)
        new = f"{uid}.pdb"
        src = os.path.join(PDB, name)
        dst = os.path.join(PDB, new)
        if os.path.exists(dst):
            # уже есть файл с правильным именем — оставляем как есть
            skipped += 1
            continue
        os.rename(src, dst)
        renamed += 1

    print(f"Переименовано: {renamed}")
    print(f"Пропущено (не подошло под шаблон или уже было): {skipped}")
    print(
        f"Всего .pdb сейчас: {len([f for f in os.listdir(PDB) if f.endswith('.pdb')])}"
    )


if __name__ == "__main__":
    main()
