# Benchmark-H1: первый воспроизводимый development pilot

4 ноября 2026. Работа изолирована в ветке benchmark/h1, отдельный worktree
EnzymeKinetics-Prediction-benchmark-h1 от main 181c45e. PR C1 №2 не изменялся,
не объединялся; Harness, архитектурные черновики и ADR-0002 сохранены.
Новых commit/push/merge нет. Физические алгоритмы и окружение не менялись.

## Доступные источники и модели

Снимок версий/файлов: [source-inventory.json](source-inventory.json). Размеры
файлов — decimal MB/GB, оценки памяти ниже отдельно; веса не скачивались.

| Источник | Проверенная версия / лицензия | Доступность / размер / ресурсы | Решение H1 |
|---|---|---|---|
| CatPred-DB | upstream 931583821ff875e873002730fbbdf1047acfde77; MIT для кода; условия исходных BRENDA/SABIO данных отдельно не подтверждены | Локальные derived trainval: kcat 20 835 строк, Km 37 056; есть исходные test и cluster subsets. Локальная копия без .git, её связь с upstream SHA не установлена. CSV по 31,5/44,0 MB, CPU | Только локальный development pilot; официальный независимый benchmark заблокирован до проверки происхождения/units/splits |
| KinBench | описан в EnzCast preprint 2026.04.28.721430v1; версия release и лицензия не подтверждены | Указанный авторами Mxc666/EnzCast возвращает 404 / repository not found. Данные/split hashes/ресурсы независимо не проверены | BLOCKED; не подменять другим KinBench/kinase benchmark |
| EnzGFM | code HEAD 0631099828e1b2aadbfe56751eec5f67a7bfdf0c; лицензия кода не подтверждена, root LICENSE=404; Zenodo 22042585 указывает CC-BY-4.0 | Архив 650M 7 162 193 759 bytes; 1.5B 5 778 657 794 bytes. README упоминает 150M, но в этом release такого архива нет. Только FP32 параметры требуют примерно 2,6/6 GB плюс runtime/активации; фактический runtime не измерен | Отложен; не загружался и не объявляется совместимым с текущим окружением |
| ESM-2 t6 8M | facebook/esm2_t6_8M_UR50D @ c731040fcd8d73dceaa04b0a8e6329b345b0f5df; MIT | model.safetensors 31 384 292 bytes, 6 слоёв/320-dim; локальный cache полный. CPU batch 2, без градиентов и сети | B1 выполнен; весовой SHA и tokenizer/config hashes сохранены |
| ChemBERTa-77M-MLM | DeepChem checkpoint ed8a5374f2024ec8da53760af91a33fb8f6a15ff; model card/license metadata отсутствуют | pytorch_model.bin 13 734 116 bytes, config 3 слоя/384-dim. 77M относится к размеру pretraining corpus, не числу параметров. Локальных весов нет; CPU frozen inference возможен по архитектуре, бюджет пока не измерен | B2 BLOCKED до уточнения лицензии конкретного checkpoint, pooling/tokenization/безопасной загрузки |
| RDKit/Morgan | установлен RDKit 2026.9.1; BSD-3-Clause | Уже используется FeatureBuilder: radius=2, 256 bits, шесть дескрипторов и invalid flag; CPU, внешних весов нет | B0/B1 выполнены |

Первичные источники: [CatPred-DB](https://github.com/maranasgroup/CatPred-DB),
[код split](https://github.com/maranasgroup/CatPred-DB/blob/931583821ff875e873002730fbbdf1047acfde77/scripts/make_splits.py),
[EnzCast/KinBench — публикация](https://www.biorxiv.org/content/10.64898/2026.04.28.721430v1.full),
[EnzGFM](https://github.com/DeepBxM/EnzGFM),
[Zenodo release metadata](https://zenodo.org/api/records/22042585),
[ESM-2 checkpoint](https://huggingface.co/facebook/esm2_t6_8M_UR50D/tree/c731040fcd8d73dceaa04b0a8e6329b345b0f5df),
[ChemBERTa checkpoint](https://huggingface.co/DeepChem/ChemBERTa-77M-MLM/tree/ed8a5374f2024ec8da53760af91a33fb8f6a15ff),
[ChemBERTa-2 paper](https://arxiv.org/abs/2209.01712),
[RDKit license](https://github.com/rdkit/rdkit/blob/master/license.txt).
MIT репозитория не использован как доказательство лицензии всех исходных данных;
лицензия DeepChem library не перенесена автоматически на конкретные веса.

## Протокол и выполненный эксперимент

Полная предварительная спецификация — [PROTOCOL.md](PROTOCOL.md),
конфигурация — [benchmark-h1-pilot.json](../../../enzyme-kinetics/configs/benchmark-h1-pilot.json).
Из каждого trainval выбраны 192 различные последовательности длиной <=256 по
SHA-256, без отбора по y. Это ограниченный, нерепрезентативный development cohort.
Официальный test не открывался и не модифицировался. Новые split не называются
официальными test или независимым подтверждением опубликованных результатов.

Существующий training.run обучает mean/Ridge/MLP, выбирая head/checkpoint по
validation. Новые development-only и split-seed отделяют model seed от split.
Seed 11/33/71, split seed 42, максимум 20 эпох; общий train-only preprocessing.
B0: текущий handcrafted FeatureBuilder без Progres и доступной геометрии.
B1 отличается добавлением frozen ESM-2. Исходные публикационные baseline v2
имеют другие данные/представления и не используются для ранжирования H1.

| Target | train / validation / calibration / development evaluation | Eval clusters |
|---|---|---|
| kcat | 85 / 34 / 36 / 37 | 36 |
| Km | 83 / 29 / 38 / 42 | 36 |

Все шесть пар partitions имеют нулевое пересечение exact sequence, supplied
40cluster и canonical enzyme/substrate pair. EC и scaffold пересекаются —
числа сохранены в [kcat audit](kcat/audit.json) и [Km audit](km/audit.json).
Нельзя утверждать EC-OOD/scaffold-OOD или независимый порог identity 40%.
Pretraining membership остаётся UNKNOWN. Неизвестные assay units/условия не
превращаются в физические defaults. Km оценивается только в dataset-native
log10 concentration units, не объявленных M.

Три заранее заданных seed объединены равными весами. Повторное обучение B0
seed=11 и reload/inference совпали побайтно по прогнозам обоих targets:
[kcat](kcat/reproducibility.json), [Km](km/reproducibility.json).

## Численные результаты

| Target / arm | log10 MAE [95% cluster CI] | RMSE | R² | Spearman | 90% interval coverage / width log10 |
|---|---|---|---|---|---|
| kcat / B0 | 1.1039 [0.8479, 1.4149] | 1.3894 | -0.4601 | -0.1253 | 97.3% / 6.0385 |
| kcat / B1 | 0.9324 [0.7098, 1.2110] | 1.2169 | -0.1199 | -0.0410 | 100.0% / 7.3380 |
| km / B0 | 0.9791 [0.7083, 1.2445] | 1.3010 | 0.0701 | 0.5604 | 69.0% / 3.2660 |
| km / B1 | 1.0026 [0.7187, 1.2729] | 1.3312 | 0.0265 | 0.3452 | 73.8% / 3.6618 |

500 bootstrap resamples по eval sequence clusters дают условные percentile
интервалы; остальные CI, каждый seed и selected head — [pilot-results.json](pilot-results.json).

- kcat: paired B1−B0 MAE = -0.1715, 95% CI [-0.2791, -0.0520].
- km: paired B1−B0 MAE = 0.0236, 95% CI [-0.0047, 0.0487].

Для kcat B1 уменьшил ошибку данного пилота, но обе R² отрицательны, а Spearman
близок к нулю. Этот результат не подтверждает переносимость или превосходство
на официальном benchmark. Для Km улучшения нет. Интервалы Km покрыли только
69–74% вместо номинальных 90%; калибровка для deployment НЕ подтверждена.
Одинаковая процедура выбора head допускает разные выбранные семейства:
сравнивается добавление embedding в адаптивный baseline, а не encoder при
фиксированном типе head. Для будущего полного ablation нужен также matched-head
контроль. CI условны на этот небольшой split, без неопределённости training-set,
pretraining corpus, assay units и provenance кластеров; это не общая проверка
статистической значимости и не основание выбирать модель по official test.

## Измеренные ресурсы и воспроизводимость

Полный run: 48,84 s; frozen ESM encoding: 11,79 s kcat + 8,25 s Km.
Пик RSS процесса: 767,6 MiB. CPU, 2 threads; GPU и сеть для кодирования не
использовались. В первый B0 seed runtime включён повторный fit для проверки
воспроизводимости. RSS — process lifetime high-water mark, не allocation каждого
head. Для масштабирования нужны отдельные measurements; большие EnzGFM архивы
и foundation fine-tuning этим бюджетом не разрешены.

Локальный полный набор находится в
`enzyme-kinetics/outputs/benchmarks/h1-pilot-v1/` этого worktree: partition.csv,
выбранные головы и preprocessing, predictions, embedding vectors, metadata и
immutable manifest. [evidence.json](evidence.json) содержит исходные и выходные
контрольные суммы; веса ESM подтверждены в sidecars
[kcat](kcat/esm2.pt.meta.json), [Km](km/esm2.pt.meta.json).
Сами веса, исходные строки и модели в предлагаемый PR не входят.

Воспроизведение требует существующего окружения и локального ESM cache;
команда и отсутствие автоматической загрузки описаны в PROTOCOL.md.
Не перезаписывать v1: повторный запуск требует нового output path.

## Проверки и следующий PR

- Ruff check и format: passed (124 Python-файла).
- py_compile четырёх изменённых/новых Python-файлов: passed.
- pytest при приёмке: **87 passed, 1 skipped, 4 deselected**; H1 cases: 8.
  Первоначальная проверка v1: 85 passed, 1 skipped, 4 deselected.
- Wheel/sdist `uv build --offline`: passed; все восемь CLI `--help`: exit 0.
- Skip: optional imageio отсутствует, окружение не переустанавливалось.
- Исключены test_native_openmm_prepared_system_in_unicode_workspace (run_md),
  test_prepared_checkpoint_continues_without_new_minimization (MD continuation),
  test_shared_md_loop_stops_on_first_failing_window (Integrator.step),
  test_embedding_converges_and_changes_energy (QM/MM SCF).
- Remote CI нового H1 не запускался: commit/push не разрешены.

Точный предлагаемый состав отдельного H1 PR — [pr-files.txt](pr-files.txt).
Содержит только additive development-only/fixed-split режим существующего
training, offline frozen ESM interface, H1 coordinator/config/tests и новые
агрегированные отчёты. Default training и старые checkpoint/import paths
сохранены; альтернативный training loop, C2 runtime, physics changes, Harness,
исторические данные, веса и outputs в PR не входят.

Рекомендация по собственным весам: небольшие target-specific heads, block
projections/fusion и calibration artifacts после проверки units/splits и
повторяемых matched ablations. Текущие pilot heads не готовы к deployment.
Не обучать сейчас foundation model или GNN; не начинать deep fine-tuning.
После отдельного H1 PR — проверить lineage/лицензии независимых данных, затем
отдельно лицензированный ChemBERTa и matched-coverage B3. Это научная ветка,
параллельная инженерной последовательности C1 → ADR-0002 → C2, без её подмены.


## Приёмка B0/B1 и PR-ready заключение

Статус: **READY_FOR_PR для ограниченного development pilot**. Это инженерная
приёмка с явно сохранёнными научными ограничениями, не разрешение deployment
или утверждение валидности полного benchmark. Состав остаётся ровно 18 путей.

После выборочного review усилена проверка identity эмбеддингов: immutable model
revision, hashes weights/config/vocab/tokenizer_config/special_tokens_map,
полный normalized sequence-set hash, pooling/normalization, algorithm ID/source
SHA и versions. H1 перед fit проверяет identity, checksum tensor payload,
sequence membership, dtype/dimension/finite values. Offline режим исключает
сетевой fallback; модель eval/frozen/no_grad. Подмена tokenizer/sequence/model/
algorithm или payload даёт отказ. Старый API чтения и checkpoints сохранены.

v1 sidecars не переписаны: они исторические и не становятся новым strict cache
задним числом. Их проверенные SHA и прежний source manifest остаются в evidence;
текущий кандидат и проверки приёмки отмечены отдельно. Производительность v1
не объявляется заново измеренной на изменённом коде. Нового обучения H1 не было;
только synthetic tests, reload/inference и tiny frozen offline cache check.

Повторно проверены все hashes v1 outputs и source CSV, partitions и привязка
192 выбранных sequence/row IDs к исходным cluster ID. Для обоих targets нет
пропусков/отрицательных ID или конфликтующих назначений exact sequence в
локальном source. **Отсутствие общих cluster ID не доказывает независимую
проверку 40% sequence identity.** Происхождение и alignment coverage не проверены.

B0/B1 и все seed используют одинаковый split.csv. У imputer/scaler общих блоков
(sequence, substrate, conditions, geometry, EC, titration) совпадают fitted
statistics/mean/scale. Матрицы целиком различаются ожидаемо из-за ESM-координат.
Preprocessing fit получает только train. Calibration не участвует в выборе head;
перед inference передаётся whitelist без y/product/target-derived columns.
Метрики ансамблей воспроизвелись до 1e-12; paired cluster-bootstrap совпал.
Небольшое расхождение CSV представления прогнозов (~6e-8) соответствует записи
float32; это не изменение модели или чисел в JSON.

### Что показывают R² и интервалы

Для kcat SST=48.9207; SSE B0=71.4287, B1=54.7879. Поэтому R² отрицателен по
определению: хуже константы, равной среднему **оценочной** выборки. Эта константа
— пояснение метрики, не доступный при обучении baseline. Причинное объяснение
слабости encoder из этих чисел не следует: cohort мал, укорочен и нерепрезентативен,
head выбирался адаптивно, ancestry/assay provenance не подтверждены.

Для Km радиусы B0/B1 1.6330/1.8309, но q90 cluster-max residual на development
равны 2.2149/2.2679 (на calibration — 1.5140/1.6022). Покрыто только 24/36 и
25/36 eval clusters. Групповой finite-sample rank проверен; ошибки индекса квантиля
не обнаружено. Наблюдается более тяжёлый хвост ошибок оценки, а не доказана
конкретная причина shift. Малое число calibration groups и непроверенная
exchangeability не позволяют обещать 90% conditional/OOD coverage.

Нельзя «исправить» это расширением интервалов до нужной доли на уже просмотренном
evaluation. Требуется заранее определить prediction unit, увеличить/аудировать
calibration groups, отделить selection от calibration, измерить coverage с
неопределённостью и проверить новый frozen protocol на новой оценочной части.
Калибровка и полноценная научная пригодность остаются **не принятыми**.

### Следующие ограниченные работы

[H1.1](PROTOCOL.md#h11--registered-next-experiment-not-executed): сначала lineage,
assay units и licenses; затем matched-head mean/Ridge, до 1024 records/target,
независимо проверенные группы, три непересекающиеся outer evaluation folds,
фиксированные seeds и отдельная calibration. Overlap training folds означает,
что их оценки не становятся независимыми репликами. Предложенный бюджет —
2 CPU threads / 8 GiB / 30 min; запуск требует отдельного задания.

[H1.2](PROTOCOL.md#h12--small-fusion-model-contract-not-implemented): frozen
ESM-2 + Morgan/RDKit + assay conditions/masks + optional validated sidecars,
малые block projections и отдельная log10 head каждого target, <=100k trainable
parameters/model. Только существующие FeatureBuilder/models/training/predict;
одинаковые данные/splits/compute budget с B0/B1 и post-hoc calibration после
фиксации predictor. Fusion-модель не реализована.

Следующее минимальное действие — ревью подготовленного H1 diff. H1.1/H1.2 не
запускать до их отдельных условий и разрешения. C1/Harness/ADR-0002/C2/main не
изменялись, никакого commit/push/merge при этой приёмке нет.
