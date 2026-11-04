import os
import pandas as pd
from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR, PROJECT_DIR


# ==========================================
# 1. ФУНКЦИЯ ОТРИСОВКИ ДЕРЕВА ПАПОК
# ==========================================
def print_tree(startpath, max_depth=2):
    ignore_dirs = {
        ".git",
        ".venv",
        "__pycache__",
        "archive",
        "data",
    }  # Датасет огромный, его не показываем в дереве
    output = []

    for root, dirs, files in os.walk(startpath):
        dirs[:] = [d for d in dirs if d not in ignore_dirs]
        level = root.replace(startpath, "").count(os.sep)

        if level >= max_depth:
            dirs[:] = []
            continue

        indent = "│   " * level + "├── "
        output.append(f"{indent}{os.path.basename(root)}/")

        subindent = "│   " * (level + 1) + "├── "
        for f in files:
            if not f.endswith((".pyc", ".log", ".gif")):
                output.append(f"{subindent}{f}")

    return "\n".join(output)


# ==========================================
# 2. ФУНКЦИЯ ВЫВОДА ДАТАСЕТА
# ==========================================
def get_dataset_sample():
    dataset_path = DATA_DIR / "CatPred-DB/data/kcat/kcat-random_trainvaltest.csv"
    if not dataset_path.exists():
        return "❌ Датасет не найден по пути: " + str(dataset_path)

    try:
        df = pd.read_csv(dataset_path, nrows=5)
        # Оставляем только самые важные колонки для читаемости
        cols_to_show = [
            "uniprot",
            "value",
            "reactant_smiles",
            "ec",
            "temperature",
            "ph",
        ]
        available_cols = [c for c in cols_to_show if c in df.columns]

        # Обрезаем длинные строки для красоты
        df_show = df[available_cols].copy()
        if "reactant_smiles" in df_show.columns:
            df_show["reactant_smiles"] = df_show["reactant_smiles"].apply(
                lambda x: str(x)[:30] + "..."
            )

        return df_show.to_markdown(index=False)
    except Exception as e:
        return f"❌ Ошибка чтения датасета: {e}"


# ==========================================
# 3. ФУНКЦИЯ ЧТЕНИЯ ФАЙЛОВ
# ==========================================
def read_file_content(filepath):
    if not os.path.exists(filepath):
        return f"// ФАЙЛ НЕ НАЙДЕН: {filepath}"
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return f.read()
    except Exception as e:
        return f"// Ошибка чтения {filepath}: {e}"


# ==========================================
# 4. ГЕНЕРАЦИЯ ОТЧЕТА
# ==========================================
def main():
    report = []

    # Заголовок
    report.append("# 📊 Отчет о состоянии проекта Enzyme-Kinetics\n")
    report.append(
        "Этот отчет сгенерирован автоматически для проверки структуры и кода.\n"
    )

    # 1. Структура
    report.append("## 📁 Структура репозитория (без папки data и .venv)\n")
    report.append("```text")
    report.append(print_tree(str(PROJECT_DIR)))
    report.append("```\n")

    # 2. Датасет
    report.append("## 📈 Первые 5 строк из датасета CatPred (kcat)\n")
    report.append(get_dataset_sample())
    report.append("\n")

    # 3. Код файлов
    report.append("## 💻 Исходный код ключевых модулей\n")

    files_to_dump = [
        (
            "src/enzymekinetics_prediction/experiments/core_chemistry.py",
            "Химический модуль (ретросинтез, EC-классы)",
        ),
        (
            "src/enzymekinetics_prediction/experiments/core_structure.py",
            "Структурный модуль (PDB, липидная капля, крышка)",
        ),
        (
            "src/enzymekinetics_prediction/experiments/core_physics.py",
            "Физический модуль (ZPE, EVB, кинетика)",
        ),
        ("apps/catalysis.py", "Главный интерфейс Streamlit"),
        (
            "scripts/demos/graph_network.py",
            "Демонстрационный forward pass без обучения",
        ),
    ]

    for filepath, description in files_to_dump:
        report.append(f"### 📄 {filepath}")
        report.append(f"**Описание:** {description}\n")
        report.append(f"```python\n{read_file_content(PROJECT_DIR / filepath)}\n```\n")

    # Сохранение в файл
    final_text = "\n".join(report)
    report_path = OUTPUT_DIR / "reports" / "project_report.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as f:
        f.write(final_text)

    print(f"✅ Отчет успешно создан: {report_path}")
    print(
        "👉 Открой его в VS Code (Markdown Preview), чтобы увидеть красивое дерево и код."
    )


if __name__ == "__main__":
    main()
