from enzymekinetics_prediction.paths import PROJECT_DIR
import os
import pandas as pd


def main():
    BASE_DIR = str(PROJECT_DIR)
    DATA_DIR = os.path.join(BASE_DIR, "data", "Baseline", "kcat")
    PDB_DIR = os.path.join(BASE_DIR, "data", "alphafold_pdb")
    os.makedirs(PDB_DIR, exist_ok=True)

    print("Читаем CSV...")
    train = pd.read_csv(os.path.join(DATA_DIR, "kcat-random_trainval.csv"))
    test = pd.read_csv(os.path.join(DATA_DIR, "kcat-random_test.csv"))

    # Собираем уникальные UniProt ID
    all_uniprots = pd.concat([train["uniprot"], test["uniprot"]]).dropna().unique()
    print(f"Всего уникальных белков: {len(all_uniprots)}")

    # Формируем список для скачивания
    missing_urls = []
    for uid in all_uniprots:
        uid = str(uid).strip()
        filepath = os.path.join(PDB_DIR, f"{uid}.pdb")

        # Если файла еще нет, добавляем в список
        if not os.path.exists(filepath):
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
            # aria2 требует формат: ссылка \n out=имя_файла
            missing_urls.append(url + "\n")
            missing_urls.append(f"  out={uid}.pdb\n")

    print(f"Нужно скачать: {len(missing_urls) // 2} файлов.")

    if missing_urls:
        with open(os.path.join(BASE_DIR, "data", "manifests", "urls.txt"), "w") as f:
            f.writelines(missing_urls)
        print("✅ Список ссылок сохранен в urls.txt")
    else:
        print("Все файлы уже скачаны!")


if __name__ == "__main__":
    main()
