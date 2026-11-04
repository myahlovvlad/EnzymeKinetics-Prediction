from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import pandas as pd


def main():
    BASE = str(PROJECT_DIR)
    DATA = os.path.join(BASE, "data", "Baseline", "kcat")
    PDB = os.path.join(BASE, "data", "alphafold_pdb")
    CIF = os.path.join(BASE, "data", "alphafold_cif")
    os.makedirs(CIF, exist_ok=True)

    train = pd.read_csv(os.path.join(DATA, "kcat-random_trainval.csv"))
    test = pd.read_csv(os.path.join(DATA, "kcat-random_test.csv"))
    uids = sorted(
        {
            str(u).strip()
            for u in pd.concat([train["uniprot"], test["uniprot"]]).dropna().unique()
        }
    )

    have_pdb = {f[:-4] for f in os.listdir(PDB) if f.endswith(".pdb")}
    have_cif = {f[:-4] for f in os.listdir(CIF) if f.endswith(".cif")}
    missing = [u for u in uids if u not in have_pdb and u not in have_cif]

    print(f"Всего уникальных UniProt ID: {len(uids)}")
    print(f"Есть PDB:                    {len(have_pdb & set(uids))}")
    print(f"Есть CIF:                    {len(have_cif & set(uids))}")
    print(f"Не хватает:                  {len(missing)}")

    path = os.path.join(BASE, "data", "manifests", "missing_uids.txt")
    with open(path, "w") as f:
        f.write("\n".join(missing))
    print(f"Список сохранён: {path}")


if __name__ == "__main__":
    main()
