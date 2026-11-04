# C1: итоговая инженерная приёмка и статистическая переоценка v2

4 ноября 2026. Ветка refactor/physics-canonicalization, HEAD/main/origin/main
181c45ec0682adcdcd7b50da745267337cf1784a. Без commit/push/merge.
Это дополнение к исходному C1-отчёту, сохранённому без изменений; его статус
«статистика требует пересчёта» относится к предыдущему проходу.

## Область и исправленный алгоритм

Одна реализация в physics/geometry.py: orthorhombic minimum image,
восстановление по полному графу белковых связей до выбора CA, ближайшие образы
целых компонент относительно первой, proper Kabsch к кадру 0. Равные веса;
RMSD/RMSF/Rg в nm; DCCM безразмерная, population normalization 1/T,
zero-variance entries = 0 по соглашению. RMSF/DCCM относительно временного
среднего после той же общей подгонки. Это компактная сборка, не temporal
unwrapping или анализ диффузии; геометрический Rg не mass-weighted.

Канонический MIC используют также glu43_geometry и pbc_clash_audit; численные
метрики dynamics/research/diagnostic используют общий structural_metrics.
Дополнительная проверка в этом проходе отклоняет неоднозначные half-box образы
при реконструкции. Косые, повёрнутые, вырожденные ячейки, разорванная топология
и незамкнутые периодические циклы явно отвергаются; triclinic не поддержан.
Детерминированные тесты покрывают переносы цепей и отдельных атомов на целые
векторы ячейки, сохранение геометрии и RMSD/RMSF/Rg/DCCM, DCD round-trip,
proper rotation, независимый аналитический reference и отказы.

Формулы поля, статистики, force field, параметры/seed исходных MD не менялись.
Поле enzyme-trace остаётся непериодическим screened Coulomb. Periodic System
требует explicit prepared_nonperiodic после внешней подготовки; служебная
ячейка DCD непериодического System не вызывает ложного отказа. Старые CLI
сохранены. Сравнение ensemble_observables с функцией main на 100 seeded
fixtures: max |difference| = 6.661338147750939e-16, единицы идентичны;
[main-equivalence.json](main-equivalence.json).

## Переоценка имеющихся данных

Все три DCD/topology/state доступны локально. 50 кадров на реплику; новых MD,
QM/MM, обучения или изменения эксперимента не было. Новый отдельный набор:
enzyme-kinetics/outputs/experiments/active-site-v3/c1-acceptance-v2/.
Исходный c1-pbc output, исторические JSON/CSV и первоначальный C1 evidence
сохранены. SHA-256 источников, кода, версий и артефактов — [evidence.json](evidence.json).
Локальный manifest заморожен до отдельного повторного запуска diagnostic-101;
его собственный SHA и дополнительные portable артефакты указаны в evidence.

| Реплика | RMSD первой / второй половины, nm | tau, кадров | ESS RMSD |
|---|---|---|---|
| 101 | 0.053592305 / 0.070747291 | 10.381712 | 4.816161 |
| 202 | 0.056961932 / 0.066864907 | 7.420538 | 6.738056 |
| 303 | 0.054998312 / 0.072951780 | 10.645045 | 4.697021 |

101: raw-image control второй половины 1.232758893 nm; старое значение
1.232758915 nm отличается округлением. Реконструированное значение совпало
с независимым diagnostic-101 до 1e-12 nm; 10 кадров смены образа.
[diagnostic-101.json](diagnostic-101.json).

RMSF, Rg, DCCM и sidechain torsions пересчитаны для каждой реплики. Здесь
хранятся residue-rmsf.csv, structural-series.csv, site-dccm.csv и sampling.json
с полной ACF, tau и ESS. Полные 250x250 DCCM и raw-image controls сохранены
локально в NPZ, их hashes включены в manifest. ACF — существующая biased
оценка, tau — существующая positive-sequence оценка, ESS=N/tau; её формула
не заменялась. include_acf только открывает уже вычисленный массив.

Два H3 сравнения: единица анализа — независимая реплика, не кадр или остаток.
Exact two-sided sign-flip и BH для двух проверок:

| Контраст RMSF | Средняя разность, nm | p | q BH |
|---|---|---|---|
| catalytic − scaffold | -0.003276450 | 0.25 | 0.50 |
| positioning − scaffold | -0.000727283 | 0.75 | 0.75 |

[summary.json](summary.json) содержит все три разности для каждого теста.
Ни один результат не значим при .05; при n=3 минимальный двусторонний p=.25.
Общая стартовая структура, короткая выборка и предположения независимости/
обменности знаков ограничивают выводы. ESS — exploratory diagnostic, не число
независимых реплик. Не подтверждены равновесие, отсутствие эффекта, каталитический
механизм, устойчивость ES или причинная связь по DCCM. Старые научные выводы
не объявляются валидными. Исторические H3 JSON и report.py export не заменены.

Пересозданы [RMSD](rmsd.png), [RMSF/температура](replicas.png),
[ACF/ESS](sampling.png) и raw/corrected DCCM для
[101](replica-101-dccm.png), [202](replica-202-dccm.png),
[303](replica-303-dccm.png). Оси — saved frame/recorded step, без догадки
о физическом интервале DCD.

## Воспроизведение

Из корня репозитория, с имеющимся uv окружением и всеми sources из evidence:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 uv run --no-sync python \
  enzyme-kinetics/scripts/diagnostics/reassess_c1_replicas.py \
  --replica enzyme-kinetics/outputs/experiments/active-site-v3/replica-101 \
  --replica enzyme-kinetics/outputs/experiments/active-site-v3/replica-202 \
  --replica enzyme-kinetics/outputs/experiments/active-site-v3/replica-303 \
  --output enzyme-kinetics/outputs/experiments/active-site-v3/c1-acceptance-v3
```

Нужен новый, ранее не существовавший output path: overwrite запрещён.
Здесь v3 — пример следующего запуска, он не выполнялся. Отсутствие любого DCD,
PDB/state.csv — BLOCKED с путём отсутствующего входа. DCD не входят в Git;
checkout с одними JSON/CSV не воспроизводит реконструкцию. В текущем окружении
пропущенных реплик нет. Перед replay сверить hashes исходников; byte-identical
PNG/NPZ между версиями библиотек не обещается, сравнивать численные массивы.

## Изолированные проверки и состав PR

C1 наложен по точному списку на отдельный detached worktree main 181c45e,
без Harness и архитектурных документов. Проверено, что импорт geometry идёт
из этого worktree; использовано существующее окружение без переустановки.
Полные команды/выходы: [validation.json](validation.json).

- Ruff check: passed; format --check: 127 files formatted.
- py_compile: все 12 Python-файлов C1 passed.
- pytest: **96 passed, 1 skipped, 4 deselected**. Skip: optional imageio отсутствует.
- uv pip check: 92 packages compatible. uv build --offline: wheel/sdist built.
- Все восемь зарегистрированных CLI: --help exit 0. CLI defaults не менялись.
- История/Harness: [preservation.json](preservation.json), 159 SHA-256 совпали.
  Исходный index пуст; HEAD не изменился; divergence с origin/main 0/0.

Четыре исключённых теста (запрет на MD/QM/MM имеет приоритет):

1. test_physics.py::test_native_openmm_prepared_system_in_unicode_workspace — run_md.
2. test_research.py::test_prepared_checkpoint_continues_without_new_minimization — run_md/checkpoint.
3. test_lowT_quench.py::test_shared_md_loop_stops_on_first_failing_window — Integrator.step.
4. test_qmmm_embedding_smoke.py::test_embedding_converges_and_changes_energy — QM/MM SCF.

Не выполнены remote CI именно этого diff и перечисленные native-тесты;
локальная проверка не закрывает полный CI gate. Workflow не менялся.
Никаких новых установок, изменений uv.lock, ADR-0001 или глобальной конфигурации.

Точный allowlist — [pr-files.txt](pr-files.txt). В него входят только C1 и его
новые результаты/минимальная документация; [change-classification.json](change-classification.json)
показывает исключённые Harness/Architecture. Other отсутствует. Review patch
создан в ignored outputs/c1-acceptance-v2/review/C1.patch; это локальный артефакт,
не commit. Его SHA и точные hashes выбранных файлов — review/manifest.json.
Готово для ревью/явного разрешения фиксации с указанными ограничениями;
полная CI/научная приёмка не объявлена. Архитектурный и C2 code в PR не входят.
