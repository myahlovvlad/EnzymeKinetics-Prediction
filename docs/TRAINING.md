# Обучение и предсказания

## Подготовка

```sh
uv sync --locked
uv run --locked enzyme-train --help
```

Ожидаемые входы для каждого target `kcat`, `km` (также поддерживается `ki`):

```text
enzyme-kinetics/data/Baseline/<target>/<target>-random_trainval.csv
enzyme-kinetics/data/Baseline/<target>/<target>-random_test.csv
enzyme-kinetics/data/CatPred-DB/data/catpred_progres_embeds_dict.pt
enzyme-kinetics/data/alphafold_pdb/<UniProt>.pdb
```

CSV обязательны. Отсутствующие эмбеддинги или структуры отмечаются в coverage;
при наличии файла эмбеддингов неподдерживаемый/небезопасный формат не обходится
автоматической небезопасной загрузкой pickle. Пустые/нечисловые targets и пустые
последовательности исключаются с записью числа строк.

Крупные данные не включены в Git. Их локальная копия существует, но её upstream
commit и лицензии всех составляющих ещё требуют проверки перед распространением.

## Команды

```sh
uv run --locked enzyme-train --target kcat --seed 42 --threads 4 --epochs 120 \
  --output enzyme-kinetics/outputs/training/kcat-new

uv run --locked enzyme-train --target km --seed 42 --threads 4 --epochs 120 \
  --output enzyme-kinetics/outputs/training/km-new
```

Каталог запуска должен быть пустым или отсутствовать: существующие результаты
не перезаписываются. Без `--output` создаётся каталог с UTC timestamp.
`--data-dir` позволяет указать другую копию данных.
Подмножество моделей: `--models mean ridge pls hist mlp radial`.
`--patience` ограничивает число эпох без улучшения validation MAE.

Максимум 120 эпох не означает, что все они выполняются: early stopping
возвращает checkpoint с лучшим validation MAE. AdamW использует lr=0.001,
weight decay=0.001; target стандартизуется по train. Нейросети обучаются с MSE,
выбираются по MAE. Количество потоков ограничено; использование GPU не требуется.

Скрипт `scripts/training/train_baseline_pls.py` сохранён как короткий совместимый
вход для mean/Ridge/PLS; старая реализация выбора лучшей модели по test удалена.

## Выбор и оценка

20% групп точных последовательностей trainval выделяются под validation.
Поэтому пропорции строк могут отличаться от 80/20. Scaler, imputer и EC vocabulary
не видят validation/test при fit. Ridge alpha и число PLS-компонент выбираются
по validation; для нейросетей так же выбирается checkpoint.

После выбора модели фиксированные варианты оцениваются на исходном test,
на test без канонических совпадений белок–субстрат с trainval и на подмножестве
точно новых последовательностей. Последнее не является 40%-cluster OOD:
родственные белки могут оставаться в обучении.

`metadata.json` фиксирует конфигурацию, версии и SHA-256 CSV/эмбеддингов;
`split.csv` — конкретные записи и группы. Новые запуски также сохраняют хеши
исходников learning-модуля в начале. Первые v1-эксперименты выполнены до
добавления этого поля и до исходного Git commit; их разбиение, веса и данные
сохранены, но точный исходный commit отсутствует.

Первичные v1-результаты — один seed, без оценки статистической значимости.
Расширение исследования: несколько seed, заранее выбранный cluster split,
ablation, внешняя проверка и калибровка неопределённости.

## Inference

```sh
uv run --locked enzyme-predict \
  --run enzyme-kinetics/outputs/training/kcat-v1 \
  --input enzyme-kinetics/data/Baseline/kcat/kcat-random_test.csv \
  --output enzyme-kinetics/outputs/predictions/kcat.csv
```

По умолчанию используется модель, выбранная по validation, даже если другая
оказалась немного лучше на test. `--model radial` позволяет явно проверить
сохранённую адаптацию исходной сети.

Вход содержит последовательность, субстрат и EC; UniProt нужен для сопоставления
имеющихся эмбеддингов и PDB. Выход содержит target, `prediction_log10`,
`prediction_dataset_units` и флаги невалидного SMILES/отсутствующих признаков.
Для kcat единица — s⁻¹; для Km вывод идёт в исходной единице датасета,
которую нельзя без сверки assay переименовывать в M. Это отдельная задача
стандартизации исходных экспериментальных данных.

Модель не оценивает достоверный барьер реакции, механизм или калиброванную
вероятность. Автоматически вычислять kcat/Km из несопоставленных assay нельзя.


## Версия 0.2: физические блоки и строгие split

```sh
uv run --locked enzyme-train --target kcat \
  --models mean ridge mbpls opls ec_mbpls mlp \
  --split-scheme grouped-three-way --group-by ec --ec-level 2 \
  --physical-table enzyme-kinetics/data/physical/pockets.csv \
  --physical-table enzyme-kinetics/data/physical/quantum.csv \
  --output enzyme-kinetics/outputs/training/my-ec-run
```

`--group-by sequence_40cluster` или `sequence_60cluster` использует одноимённую
колонку существующего датасета; отсутствующие группы отвергаются.
`grouped-three-way` создаёт новый test по группам из объединённых исходных CSV,
затем независимую validation по оставшимся группам. Это другой эксперимент,
а не прежняя official test-таблица. `split.csv` содержит partition,
source_dataset/source_record и split_group; metadata — число пересечений.
Происхождение и качество сторонних cluster assignments требуют отдельного аудита.

Физические таблицы требуют sidecar/checksum и не содержат targets.
[Формат и генерация](PHYSICS.md). `--physical-table` повторяется для каждого
блока. Все модели получают одинаковые входы; локальные EC-проекции используют
только train. Новая EC-группа автоматически получает глобальный прогноз.

Для настоящих ESM-2-векторов:

```sh
uv run --locked --extra embeddings enzyme-embed \
  --input enzyme-kinetics/outputs/backend_checks/esm-input.csv \
  --output enzyme-kinetics/outputs/embeddings/my-esm2.pt
uv run --locked enzyme-train --target kcat \
  --embedding-path enzyme-kinetics/outputs/embeddings/my-esm2.pt
```

Вход CSV должен содержать sequence; примерный путь — локальный двухстрочный
smoke input, для полного обучения нужен полный набор последовательностей.
Exporter фиксирует model revision, residue-mean pooling и SHA-256, сохраняет
ключи exact sequence; слишком длинные последовательности отвергаются,
тихого усечения нет. Двухстрочный пилот не заменяет покрытие всего датасета.
Сравнения v2 и ограничения: [отчёт](benchmarks/v2/README.md).


## Версия признаков 3: pI и титрование

Новые FeatureBuilder.fit добавляют блок `titration`: sequence pI как корень
Q(pH), модельный заряд при pH assay и pH−pI. Используется модель Biopython;
неизвестные аминокислоты/условия дают явные NaN, которые обрабатываются
train-only imputer. Это не измеренный pI и не расчёт каталитических pKa.
Старые preprocessing v1/v2 сохраняют прежнюю размерность через версионный
атрибут builder. Обучение kcat-titration-validated-v3 выполнено; результаты
и проверки загрузки находятся в docs/experiments/active-site-v3/.
