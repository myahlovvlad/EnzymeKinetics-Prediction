# C1: PBC и структурные метрики

Дата протокола: 4 ноября 2026. База: 181c45e; незакоммиченный C1 в
refactor/physics-canonicalization. Это повторный анализ готовых траекторий,
не новый эксперимент MD и не доказательство устойчивости ансамбля.

## Контракт и проверяемая гипотеза

Гипотеза: скачок RMSD replica 101 возникает из-за смены периодического образа
цепи B. Связанная реконструкция перед общей CA-подгонкой должна воспроизвести
ранее диагностированные 0.05359/0.07075 nm вместо 0.05359/1.23276 nm.
Критерий: полусредние совпадают с прежней диагностикой с абсолютным допуском
1e-5 nm (исторический MDTraj RMSD float32 против текущего Kabsch float64).
Это не критерий научной приёмки MD или кинетической модели.

Каноническая реализация — physics/geometry.py, метод
orthorhombic-bonded-assembly-kabsch-v1 (владение physics согласно ADR-0001).
Условия:

- Координаты/длины/векторы ячейки — nm; время DCD — ps. Векторы — строки.
- Поддерживаются положительные прямоугольные ячейки, ориентированные по осям,
  со своей матрицей для каждого кадра. Косые/повёрнутые/вырожденные ячейки
  явно отвергаются; одни длины без углов не используются вместо матрицы.
- Полный белковый граф связей реконструируется до CA-выборки. Связанные
  компоненты делаются целыми; незамкнутые периодические циклы и неполная
  связность цепи отвергаются. Связи не угадываются по расстояниям.
- Целые молекулы помещаются ближайшими геометрическими центрами к первой
  связанной компоненте. Это допущение компактной сборки из прежней диагностики
  димера, не temporal unwrapping и не анализ диффузии отдельных молекул.
- Все выбранные CA подгоняются с одинаковыми весами к кадру 0: перенос и
  proper Kabsch rotation без отражения. RMSD относительно кадра 0, RMSF —
  относительно временного среднего после этой же подгонки, Rg геометрический.
  DCCM — средние скалярные произведения флуктуаций / sqrt дисперсий, делитель T.
  При нулевой дисперсии DCCM=0 по явному соглашению (не измеренная корреляция).
- Структурные массивы отделены от autocorrelation/ESS, температурной статистики,
  ролей остатков и межрепличных проверок. Масса атомов не входит в эти метрики.
- Непериодическое screened-Coulomb поле не заменяется PME/MIC-суммой. Для
  DCD с ячейкой enzyme-trace требует явного --prepared-nonperiodic после
  отдельной подготовки образов всех нужных источников и probe-атомов.

## Источники и воспроизведение

В этом окружении найдены все три пары trajectory.dcd/topology.pdb, а также
state.csv, в enzyme-kinetics/outputs/experiments/active-site-v3/replica-{101,202,303}.
У каждой траектории 50 кадров, 46690 атомов, 3760 белковых атомов и две
связанные белковые компоненты по 1880 атомов. В Git хранятся исторические
отчёты, но не сами DCD: в другом checkout отсутствие источников — BLOCKED,
а не повод считать старый JSON пересчитанным.

Параметры моделирования/force field не меняются. Seed анализа отсутствует:
алгоритм детерминирован; synthetic regression использует 101/202.
Новые локальные результаты записываются в outputs/experiments/active-site-v3/c1-pbc/
с отдельными replica-101/202/303 и diagnostic-101.json. Непустые выходы
отвергаются. Прежние JSON/CSV/DCD не перезаписываются.

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 uv run --no-sync python - <<'PY'
from pathlib import Path
from enzymekinetics_prediction.research.ensemble import analyze_replica
root = Path('enzyme-kinetics/outputs/experiments/active-site-v3')
for seed in (101, 202, 303):
    analyze_replica(root / f'replica-{seed}', root / 'c1-pbc' / f'replica-{seed}')
PY
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 uv run --no-sync python \
  enzyme-kinetics/scripts/diagnose_dimer_rmsd.py \
  --trajectory enzyme-kinetics/outputs/experiments/active-site-v3/replica-101/trajectory.dcd \
  --topology enzyme-kinetics/outputs/experiments/active-site-v3/replica-101/topology.pdb \
  --output enzyme-kinetics/outputs/experiments/active-site-v3/c1-pbc/diagnostic-101.json
```

В новых analysis.json фиксируются метод, версии NumPy/MDTraj, хеши исходных
DCD/топологии и трёх модулей анализа; structure.npz хранит RMSD/RMSF/Rg/DCCM.
Компактное portable evidence с хешами и результатами — [evidence.json](evidence.json).
Старый H3-replicated-md.json и его p/q не обновляются:
они требуют отдельного согласованного межрепличного статистического анализа.


## Выполненный повторный анализ

| Реплика | Исторический RMSD: первая / вторая половина, nm | C1: первая / вторая половина, nm |
| --- | --- | --- |
| 101 | 0.053592308 / 1.232758915 | 0.053592305 / 0.070747291 |
| 202 | 0.056961930 / 0.066864904 | 0.056961932 / 0.066864907 |
| 303 | 0.054998318 / 0.072951770 | 0.054998312 / 0.072951780 |

101: 10 кадров со сменой периодического образа. Ранее опубликованная
диагностика воспроизведена в пределах 1e-5 nm; новый research и диагностический
пути совпадают до 1e-12 nm. Для 202/303 расхождения на уровне округления.
Новые RMSF и DCCM также рассчитаны и сохранены в отдельных CSV/structure.npz;
старые файлы не исправлялись на месте. Например, для 101 контраст
positioning-minus-scaffold RMSF меняется с 0.110073 до -0.000657 nm, а
catalytic-minus-scaffold — с 0.009531 до -0.001948 nm. Это описательные
пересчитанные значения; знак старого эффекта не переносится в новый вывод.

**Не подтверждено:** равновесие/сходимость MD, валидность ES или механизма,
причинность корреляций. Старые H3 p/q, объединённые графики и интерпретации,
построенные по прежнему RMSF/DCCM 101, требуют отдельного согласованного
пересчёта статистики/визуализации всей серии. Оценки ESS остаются лишь
exploratory short-series diagnostics; дополнительные независимые реплики
и длинное семплирование в C1 не выполнялись.

**Наличие данных:** все DCD/topology/state для 101/202/303 доступны локально,
пропущенных реплик здесь нет. В Git DCD не поставляются. При отсутствии
любого источника нужно отметить BLOCKED и проверить checksum; CSV/JSON
истории недостаточно для повторной реконструкции. В evidence перечислены
ожидаемые пути и SHA-256, включая state.csv и код статистики/аннотаций.

## Проверки и границы доказательств

- 12 новых детерминированных тестовых случаев: целые векторы ячейки для цепей
  и отдельных атомов, сохранение геометрии, реконструкция по связям до CA,
  RMSD/RMSF/DCCM/Rg, аналитический tetrahedron reference, proper rotation,
  отказ при косой ячейке/неполной топологии/незамкнутом цикле, синтетический
  DCD round-trip, сохранение старых результатов, отказ неподготовленного trace.
- Pytest: **96 passed, 1 skipped, 4 deselected**. Skip — optional imageio.
  Исключены два run_md/checkpoint теста, тест с Integrator.step и QM/MM
  embedding smoke: они нарушили бы явный запрет этой задачи на MD/QM/MM.
- Ruff check/format, py_compile восьми C1 Python-файлов, uv pip check и
  uv build --offline (wheel/sdist) успешны. Зависимости не переустанавливались.
- SHA-256 112 ранее существовавших файлов совпали, включая исторические
  JSON/CSV, исходные DCD/топологии и незакоммиченные файлы harness.
- Workflow GitHub Actions не менялся. CI именно этого незакоммиченного C1
  **не запускался**; предыдущий зелёный main не доказывает успех нового CI.
  Полный набор native-MD/QM/MM тестов этим запуском не проверен.

Воспроизведение разрешённой проверки без моделирования:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 uv run --no-sync pytest -q -ra \
  --deselect enzyme-kinetics/tests/test_physics.py::test_native_openmm_prepared_system_in_unicode_workspace \
  --deselect enzyme-kinetics/tests/test_research.py::test_prepared_checkpoint_continues_without_new_minimization \
  --deselect enzyme-kinetics/tests/test_lowT_quench.py::test_shared_md_loop_stops_on_first_failing_window \
  --deselect enzyme-kinetics/tests/test_qmmm_embedding_smoke.py::test_embedding_converges_and_changes_energy
uv run --no-sync ruff check enzyme-kinetics
uv run --no-sync ruff format --check enzyme-kinetics
uv pip check
uv build --offline
```

C1 меняет только подготовку координат/анализ и защиту новых выходов; параметры
силового поля, интегратор и физические формулы электрического поля сохранены.
Повторный вызов analyze_replica для исторического run пишет в новый analysis-c1,
а не заменяет его analysis.json. Существующий экспорт research/report.py всё
ещё читает исторический корневой analysis.json: автоматическое смешение старой
статистики с новыми CSV намеренно не выполняется. Новый отчёт нужно строить
явно по согласованному комплекту C1-результатов.
