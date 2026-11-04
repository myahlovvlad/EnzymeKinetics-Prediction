# Физический вычислительный контур 0.2

Контур выполняет реальные расчёты и хранит их отдельно от иллюстраций
`experiments/` и `apps/`. Результаты — исследовательские признаки с единицами,
исходными файлами, настройками, версиями, SHA-256 и ограничениями.
Работающий backend не доказывает физическую сходимость результата.

## Установка и первый запуск

```sh
uv sync --locked --extra physics --extra embeddings
uv run --locked --extra physics enzyme-physics \
  --config enzyme-kinetics/configs/p00947-apbs-check.yaml \
  --output enzyme-kinetics/outputs/physical/my-p00947-run
```

Потребуются локальная структура P00947 и распакованный P2Rank в `vendor/p2rank`.
Готовые данные проекта не входят в wheel. Проверенный P2Rank 2.4.2 использует
Java 21 из `jdk4py`; Java 25 с его Groovy несовместима. Явный `JAVA_HOME`
оператора имеет приоритет. APBS — отдельная внешняя программа. Официальная сборка 3.5.0 установлена
локально в vendor с проверкой опубликованной SHA-256. Воспроизводимая установка:

```sh
uv run --locked python enzyme-kinetics/scripts/structure/install_apbs.py
```

Приоритет поиска: `APBS_EXECUTABLE` → системный PATH → закреплённый локальный
`vendor/apbs/APBS-3.5.0.Linux/bin/apbs`. Сборка рассчитана на Linux x86_64.
На других платформах задайте собственный APBS executable.
Если выбран APBS, отсутствие executable завершает соответствующую стадию
ошибкой, а не заменяет метод автоматически.

Выходной каталог должен быть новым или пустым. `manifest.json` обновляется
после каждой стадии; ошибка записывается вместе с типом и сообщением. CLI
возвращает ненулевой код, если хотя бы одна включённая стадия не выполнена.
Успешные независимые стадии сохраняются. Подпроцессы ограничены timeout.

Для явно выбранного screened-Coulomb варианта сохранён
`configs/p00947-backend-check.yaml`; APBS-вариант — `p00947-apbs-check.yaml`.

Пилотный YAML задаёт **предполагаемую** ионную силу 0.15 M и отключает расчёт
скорости: для P00947 нет проверенного ΔG‡. Время production MD составляет
0.2 ps, пять кадров. Это проверка установки, записи DCD и извлечения признаков;
система не успевает даже выйти на заданную температуру. Это не равновесная MD,
не ES-комплекс и не данные для научного вывода о динамике катализа.

## Методы и единицы

| Блок | Реализация | Ограничение |
| --- | --- | --- |
| Карман | P2Rank, top rank, состав остатков, Rg, центр | Карман связывания не автоматически каталитический; charge proxy не электростатика |
| QM | tblite GFN2-xTB на ETKDGv3/MMFF-конформере; энергия, орбитали, заряды, диполь | Single point, не оптимизация TS, не ΔG‡; сольватация только при явной настройке |
| Fukui | Вертикальные N±1 состояния при одной геометрии | Минимальный спин по чётности; нельзя подменять исследование электронных состояний |
| Диссоциация | Разность электронных энергий явно заданных фрагментов | Не термохимическая BDE, не свободная энергия активации |
| PQR | PDB2PQR AMBER + PROPKA при заданном pH | Нужна проверка каталитических протонированных форм; нулевые радиусы полярных H допустимы в AMBER |
| APBS | LPBE input, subprocess, OpenDX potential → поле | Реальный LPBE-прогон выполнен; нужна сходимость сетки и проверка диэлектрической модели |
| Приближённое поле | Явно выбранный Debye–Hückel screened Coulomb | Однородный диэлектрик; это не решение PB с границей белок–раствор |
| MD | OpenMM CPU, Amber14/TIP3P-FB или подготовленные System XML + PDB | Стандартный apo-белок; для лигандов и металлов нужна параметризация |
| Динамические признаки | Выравнивание CA, RMSD/Rg/DCCM, заданные расстояния | Короткие моменты не устанавливают сходимость и причинную роль в катализе |
| Кинетика | Eyring, Wigner/Bell-parabolic, явный микрокинетический цикл | Только обоснованные входные барьеры и элементарные константы |

Входы QM: Å, заряд в e; внутри tblite — атомные единицы.
Выходы: Hartree/eV/Debye. Поле: V/Å; потенциал: V.
APBS/OpenDX kBT/e переводится в V при явно заданной температуре.
OpenMM/траектории: nm, ps, K, mol/L. Обучающие таблицы используют температуру
в °C, `ph`, `ionic_strength_M`; преобразование Kelvin → Celsius явное.
Km модели остаётся в единицах исходного датасета, пока assay не сверены.

Разъединённые SMILES не объявляются единой QM-геометрией. Их необходимо
рассчитывать как отдельные частицы либо задавать физический комплекс.
Неудачные расчёты сохраняются в списке failures; отсутствующие значения —
NaN/None с индикаторами пропусков, а не нулевые физические величины.

## Траектория атомов и поле вдоль выбранной пары

```sh
uv run --locked --extra physics enzyme-trace \
  --trajectory enzyme-kinetics/outputs/physical/my-p00947-run/md/trajectory.dcd \
  --topology enzyme-kinetics/outputs/physical/my-p00947-run/md/topology.pdb \
  --system enzyme-kinetics/outputs/physical/my-p00947-run/md/system.xml \
  --pair 4 21 --source-selection protein \
  --label 'Проверка CA-пары, не каталитическая связь' \
  --ionic-strength-M 0.15 --dielectric 78.5 \
  --output enzyme-kinetics/outputs/physical/my-p00947-trace
```

Индексы — нумерация **всех атомов topology с нуля**. Пример использует CA
первых двух остатков подготовленного пилота; для другой topology индексы
нужно определить заново. Для реакции указываются химически обоснованные атомы
субстрата/катализатора в параметризованном комплексе. Оба probe-атома исключаются
из источников поля, потенциал берётся в середине пары, проекция положительна
в направлении первого атома ко второму. `trace.csv` содержит каждый кадр,
расстояние и поле. Заряды извлекаются из соответствующего OpenMM System;
порядок атомов должен совпадать с trajectory/topology.

Поле анализа — непериодическое однородное screened Coulomb; оно отличается от
PME, которым могла генерироваться траектория. Минимальные периодические образы
не суммируются. Для пересечений границ нужна отдельная подготовка траектории.
Нельзя трактовать этот анализ как квантовую реакцию или APBS по кадрам.

## Протекание обратимой реакции

```sh
uv run --locked enzyme-cycle \
  --config enzyme-kinetics/configs/reversible-cycle-example.json \
  --output enzyme-kinetics/outputs/reactions/my-numerical-check
```

Это явно обозначенный **синтетический численный пример**, не параметры P00947.
Механизм: `E+S ⇌ ES ⇌ EP ⇌ E+P`. Интегратор BDF сохраняет общий фермент
`E+ES+EP` и количество субстратных молекул `S+P+ES+EP`; отрицательные
концентрации сверх допуска приводят к ошибке. Выходы: progress.csv и evidence.
Используется приближение идеального разбавленного раствора.

При `a=kon(S), b=koff(S), c=kchem_forward, d=kchem_reverse, e=release(P)`:
`kcat=c*e/(c+d+e)`, `Km=[b*(d+e)+c*e]/[a*(c+d+e)]` при начальном P≈0.
Продуктовая ассоциация учитывается на полной временной кривой.
`Keq=(kon(S)*c*release(P))/(koff(S)*d*kon(P))`.
Если задана независимая стандартная ΔG реакции и T, проверяется
`ΔG°=-RT ln(Keq)`; несогласованные константы отвергаются.
Набор kcat/Km сам по себе не идентифицирует все элементарные константы.

Туннельные поправки требуют мнимой частоты соответствующего TS и области
применимости приближения: Bell-parabolic `x/sin(x)`, Wigner `1+x²/6`,
`x=h*c*ν/(2*kBT)` и здесь `x<π`. Для KIE нужны отдельные изотопные
свободные энергии/частоты; универсальная константа по одной массе не вводится.

## Физические признаки для ML

```sh
uv run --locked --extra physics enzyme-features pockets \
  --output enzyme-kinetics/data/physical/pockets.csv
uv run --locked --extra physics enzyme-features quantum \
  --dataset enzyme-kinetics/data/Baseline/kcat/kcat-random_trainval.csv \
  --max-molecules 32 --cache-dir enzyme-kinetics/outputs/quantum-cache \
  --output enzyme-kinetics/data/physical/quantum.csv
uv run --locked enzyme-train --target kcat --models mean ridge mbpls opls ec_mbpls mlp \
  --physical-table enzyme-kinetics/data/physical/pockets.csv \
  --physical-table enzyme-kinetics/data/physical/quantum.csv \
  --output enzyme-kinetics/outputs/training/my-physical-run
```

CSV требует sidecar `.csv.meta.json`: блок, join keys, единицы каждого признака,
SHA-256, происхождение и явный флаг label-independent extraction.
Разрешены блоки pocket/qm/electro/md; feature names обязаны иметь префикс блока.
Checksum не совпадает — обучение останавливается. Lookup сохраняется вместе
с preprocessing: inference не перечитывает изменённую исходную таблицу.
Электростатика/MD связываются также по pH/T/ионной силе. Неизвестная ионная
сила assay не совпадает с предположенной 0.15 M автоматически.

## Источники методов

- [Официальный релиз APBS 3.5.0](https://github.com/Electrostatics/apbs/releases/tag/v3.5.0);
  SHA-256 Linux zip: `b92f6a1bc05c9d4278d85a68072ab6b1ce5aed8d4bb53870746104635b923ca3`.

- [tblite Python API](https://tblite.readthedocs.io/en/latest/api/python.html).
- [OpenMM: running simulations](https://docs.openmm.org/latest/userguide/application/02_running_sims.html).
- [PDB2PQR usage](https://pdb2pqr.readthedocs.io/en/latest/using/index.html).
- [APBS source/documentation](https://github.com/Electrostatics/apbs).
- [Область применимости туннельных приближений](https://pmc.ncbi.nlm.nih.gov/articles/PMC2854803/).
