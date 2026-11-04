# EnzymeKinetics-Prediction

Исследовательский проект для предсказания `kcat`/`Km` и изучения физических
признаков ферментативного катализа на CPU. Версия ПО **0.3.0** объединяет
обучение, хемометрику и исполняемый физический контур.

Реализованы реальные P2Rank, GFN2-xTB, PDB2PQR и OpenMM, анализ траекторий,
приближённое поле и обратимая микрокинетическая симуляция. APBS установлен локально из официальной сборки; расчёт Пуассона–Больцмана
также проверен. Поиск TS, QM/MM-свободные
энергии и экспериментально валидированный расчёт kcat остаются научными
задачами: HOMO–LUMO gap и короткая apo-MD их не заменяют.

**[Насколько проект соответствует идее и что делать дальше](docs/IMPLEMENTATION_STATUS.md)** ·
[План реализации](docs/specs/PHYSICAL_PIPELINE.md) · [SDLC](docs/SDLC.md)

## Активный центр и серия экспериментов 0.3

Добавлены каталитические роли и контактные остатки, chi1/chi2, pI/Q(pH),
факторный план pH×T×I, три MD-реплики WT-димера, торсионное сопряжение,
контроль сеточной сходимости и журнал с зарегистрированными гипотезами.
[Результаты и дневник](docs/experiments/active-site-v3/README.md),
[протокол](docs/experiments/PROTOCOL.md).
Полные химические TS и экспериментально проверенный kcat пока не получены.

## Быстрый старт

Один корневой uv-проект: Python 3.13, PyTorch CPU, `uv.lock`, `.venv`.

```sh
uv sync --locked
uv run --locked pytest -q
uv run --locked enzyme-train --target kcat
uv run --locked enzyme-train --target km
```

Нужны локальные CSV `enzyme-kinetics/data/Baseline/{kcat,km}/`.
Данные, vendor и обученные веса не включаются в Git.
[Входные данные](docs/DATA_LAYOUT.md), [обучение](docs/TRAINING.md).

Физический пилот с явно обозначенными допущениями:

```sh
uv sync --locked --extra physics
uv run --locked --extra physics enzyme-physics \
  --config enzyme-kinetics/configs/p00947-apbs-check.yaml \
  --output enzyme-kinetics/outputs/physical/my-p00947-run
```

Понадобятся P00947.pdb и P2Rank из локальных данных/vendor.
Пилот выполняет 0.2 ps apo-MD для проверки backend, а не достаточное
семплирование реакции. Ионная сила предположена; барьер и kcat не выдумываются.
[Протоколы, поля, траектории и кинетические кривые](docs/PHYSICS.md).

## Выполненные эксперименты

Модели выбираются по validation. Ошибки в log10; для Km исключены пять
совпадающих канонических пар, остаётся 4 113 test-строк. Один seed: 42.

| Эксперимент | Выбранная модель | MAE test | R² |
| --- | --- | --- | --- |
| kcat baseline v1 | MLP | 0.9609 | 0.3677 |
| kcat + pocket/QM v2 | MLP | 0.9883 | 0.3481 |
| Km baseline v1 | MLP | 0.6853 | 0.5031 |
| Km + pocket/QM v2 | MLP | 0.6797 | 0.5119 |
| kcat EC-OOD v2 | OPLS | 1.2451 | 0.0731 |
| kcat cluster40 v2 | MLP | 1.2057 | 0.1227 |

Физические признаки не дали устойчивого улучшения. Pocket-покрытие около
четверти строк; QM-покрытие низкое. EC-OOD — другой test с разнесёнными
подклассами EC, его числа нельзя напрямую сравнивать как улучшение v1.
Прежний test уже просматривался: новые результаты являются исследовательскими.
[Полные v2-сравнения и evidence](docs/benchmarks/v2/README.md),
[v1 kcat](docs/benchmarks/kcat/README.md), [v1 Km](docs/benchmarks/km/README.md).

В локальных `outputs/training/` сохранены модели, preprocessing, split,
предсказания и хеши. Пример inference:

```sh
uv run --locked enzyme-predict \
  --run enzyme-kinetics/outputs/training/kcat-physical-v2 \
  --input enzyme-kinetics/data/Baseline/kcat/kcat-random_test.csv \
  --output enzyme-kinetics/outputs/predictions/kcat.csv
```

## Структура

```text
pyproject.toml, uv.lock, .python-version  # единое окружение
.github/workflows/ci.yml                 # core/physics проверки
docs/                                   # спецификация, ADR, методы, evidence
archive/                                # уникальные исторические материалы
enzyme-kinetics/
  src/enzymekinetics_prediction/
    learning/                           # признаки, модели, train/inference, ESM
    physics/                            # QM, карманы, поля, MD, кинетика
    experiments/                        # прежние иллюстративные модели
    paths.py                            # общие пути
  configs/                              # проверяемые протоколы и численный пример
  scripts/                              # data, structure, diagnostics, training
  apps/                                 # иллюстративный Streamlit
  tests/                                # unit/integration/native проверки
  docx/                                 # постановка задачи и пилот P00947
  data/, outputs/, vendor/              # локальные артефакты вне Git
```

## Документация

- [Соответствие идее, границы и следующий научный цикл](docs/IMPLEMENTATION_STATUS.md)
- [Физические расчёты и команды](docs/PHYSICS.md)
- [Окружение](docs/ENVIRONMENT.md), [обучение](docs/TRAINING.md)
- [Сравнение подходов и источники](docs/APPROACHES.md)
- [Архитектура](docs/ARCHITECTURE.md), [SDLC](docs/SDLC.md), [разработка](CONTRIBUTING.md)
- [Данные](docs/DATA_LAYOUT.md), [исследовательские материалы](enzyme-kinetics/docx/README.md)
- [Реорганизация](docs/REORGANIZATION.md), [история изменений](CHANGELOG.md)
