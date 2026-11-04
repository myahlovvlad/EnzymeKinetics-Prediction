# Рабочее окружение

## Единый проект uv

Канонические файлы находятся в корне: `pyproject.toml`, `uv.lock`,
`.python-version`. Пакет собирается из `enzyme-kinetics/src/` через uv_build.
Выбран Python 3.13; проверенная локальная версия — 3.13.16.

```sh
uv sync --locked
uv run --locked python enzyme-kinetics/scripts/diagnostics/check_environment.py --profile chemistry --profile quantum --profile dashboard --profile ml
uv run --locked pytest -q
uv run --locked ruff check enzyme-kinetics
uv run --locked ruff format --check enzyme-kinetics
uv build
```

Torch, RDKit, PySCF, библиотеки анализа данных и интерфейса входят в основные
зависимости. Отдельные флаги `--extra chemistry`, `--extra quantum`,
`--extra dashboard` и `--extra ml` больше не нужны и не определены.
PyTorch разрешается через отдельный официальный CPU-index;
CUDA-пакеты не требуются. Выбор источника закреплён в манифесте.

Для ноутбуков:

```sh
uv sync --locked --extra notebooks
uv run --locked --extra notebooks jupyter lab
```

`uv sync` и `uv run` синхронизируют окружение с выбранными зависимостями;
нужные необязательные extras следует передавать повторно. См.
[официальную документацию uv](https://docs.astral.sh/uv/concepts/projects/sync/).

## IDE и пути

Открывайте корень репозитория в VS Code. Интерпретатор:
`${workspaceFolder}/.venv/bin/python`; обнаружение тестов настроено на
`enzyme-kinetics/tests`. Для импортов задан `enzyme-kinetics/src`.

Данные, результаты и vendor находятся внутри `enzyme-kinetics/`.
Общие пути задаёт `enzymekinetics_prediction.paths`; смена cwd их не меняет.
При установке wheel вне репозитория задайте переменную
`ENZYME_KINETICS_PROJECT_DIR` с путём к каталогу данных проекта.

Прежнее повреждённое окружение и промежуточное окружение Python 3.12 удалены.
На промежуточном 3.12 наблюдался аварийный импорт pandas/PyArrow;
его низкоуровневая причина не установлена. На итоговом системном Python 3.13
импорты проверены. Этот факт не является диагнозом всех установок Python 3.12.

## Дополнительные сценарии

```sh
uv run --locked python 'archive/transition_state/Работающий код .py'
uv run --locked streamlit run enzyme-kinetics/apps/catalysis.py
uv run --locked python enzyme-kinetics/scripts/data/audit_enzyme_data.py --help
uv run --locked python enzyme-kinetics/scripts/data/download_alphafold_pdbs.py --help
```

Первый скрипт создаёт молекулы RDKit/PySCF без расчёта переходного состояния.
Интерфейс требует локальную структуру Q05469 и показывает иллюстрацию.
Массовые загрузки запускаются только явным вызовом соответствующего скрипта;
импорт утилит загрузки больше не начинает скачивание.

P2Rank сохранён в `enzyme-kinetics/vendor/p2rank/`. В версии 0.2 доступны:

```sh
uv sync --locked --extra physics --extra embeddings
uv run --locked --extra physics pytest -q
```

Physics устанавливает tblite, OpenMM, MDTraj, PDB2PQR/PROPKA и Java 21
через jdk4py. Отдельный xTB executable не требуется tblite API. Эти backend
проверены реальными локальными расчётами. APBS — внешняя программа: официальная Linux-сборка 3.5.0 установлена
локально с проверкой опубликованной SHA-256; используется без изменения
системного PATH. Повторная установка: `uv run --no-sync python
enzyme-kinetics/scripts/structure/install_apbs.py`. Vina не установлена.
При отсутствии APBS адаптер выдаёт явную ошибку.
Embeddings добавляет transformers для официального ESM-2.
[Протоколы и ограничения](PHYSICS.md).

[Отчёт реорганизации](REORGANIZATION.md) содержит карту переносов.
