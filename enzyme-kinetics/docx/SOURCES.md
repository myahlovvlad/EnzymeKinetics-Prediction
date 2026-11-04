# Источники и статус проверки

Проверка выполнена 1 ноября 2026. Ниже первичные статьи, официальная документация
и базы. Это целевой обзор для постановки задачи, а не систематический обзор
всей литературы. Полные тексты статей не копировались в репозиторий.

| Источник | Что использовано | Граница проверки |
| --- | --- | --- |
| [IUBMB: Symbolism and Terminology in Enzyme Kinetics](https://iubmb.qmul.ac.uk/kinetics/ek4t6.html) | Определения наблюдаемых кинетических величин | Не задаёт механизм выбранного фермента |
| [RCSB 8CHO](https://www.rcsb.org/structure/8CHO), Cho et al., 1998, DOI 10.1021/bi9801614 | P00947, WT, димерная сборка, каталитические остатки | Карточка и аннотация; координаты не скачаны |
| [RCSB 1QJG](https://www.rcsb.org/structure/1QJG), Cho et al., 1999, DOI 10.1074/jbc.274.46.32863 | Комплекс с аналогом и наличие мутации | Точная мутация и пригодность координат требуют проверки |
| [BRENDA EC 5.3.3.1](https://brenda-enzymes.de/enzyme.php?ecno=5.3.3.1) | Совпадение Km 0.1528, условий и WT с P00947 | Кураторская запись, а не проверенная исходная таблица assay |
| [Yun et al., 2003: Origin of the Different pH Activity Profile…](https://doi.org/10.1074/jbc.M302166200) | Кандидат-первоисточник BRENDA 652424, PMID 12734184 | Библиография и описание найдены; полная таблица, буфер и связь обеих целей с assay пока не проверены |
| [Fried, Bagchi, Boxer, 2014, Science, DOI 10.1126/science.1259802](https://pubmed.ncbi.nlm.nih.gov/25525245/) | Экспериментальная мотивация поля как наблюдаемой характеристики KSI | Не подтверждает универсальность признака и не заменяет расчёт P00947 |
| [Transition state theory for enzyme kinetics](https://pmc.ncbi.nlm.nih.gov/articles/PMC4555010/) | TST, реакционная координата, ядерные эффекты и transmission | Использованы индексируемые фрагменты; прямое открытие PMC встретило проверку браузера |
| [A practical guide to modelling enzyme-catalysed reactions, 2012, DOI 10.1039/C2CS15297E](https://pubs.rsc.org/en/content/articlehtml/2012/cs/c2cs15297e) | Различие потенциальных барьеров, свободной энергии и экспериментальной скорости | Доступны индексируемые фрагменты; повторное открытие вернуло 403 |
| [PySCF: QM/MM methods](https://pyscf.org/user/qmmm.html) | Возможность embedding внешними зарядами | Не подтверждает готовность полного QM/MM workflow проекта |
| [OpenMM: Standard Forces](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html) | Формы стандартных force-field terms | Двигатель ещё не выбран и не установлен |
| [PLUMED: Using restraints](https://www.plumed.org/doc-v2.9/user-doc/html/lugano-2.html) | Umbrella sampling и снятие смещения | Ссылка на доступную версию учебного материала, не готовый KSI-протокол |
| [PyMBAR timeseries](https://github.com/choderalab/pymbar/blob/main/pymbar/timeseries.py) | Автокорреляции и эффективные выборки | Диагностика не доказывает полноту семплирования |
| [Northrup, Allison, McCammon, 1984, DOI 10.1063/1.446900](https://doi.org/10.1063/1.446900) | Броуновская ассоциация как отдельная кинетическая задача | Полный текст не изучен; используется библиографический ориентир |
| [Interfacial activation of M37 lipase, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5287222/) | Специфичность межфазного сценария | Не переносится напрямую на KSI или Q05469 |
| [Westerhuis, Kourti, MacGregor, 1998: Multiblock PCA/PLS](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291099-128X%28199809/10%2912%3A5%3C301%3A%3AAID-CEM515%3E3.0.CO%3B2-S) | Необходимость точной формулировки MB-PLS и масштабирования | Использован abstract; алгоритм для проекта ещё не выбран |
| [scikit-learn: Common pitfalls](https://scikit-learn.org/stable/common_pitfalls.html) | Предобработка без утечек | Конкретный split нужно определить под исследовательский вопрос |
| [CatPred, Boorla & Maranas, 2025, DOI 10.1038/s41467-025-57215-9](https://www.nature.com/articles/s41467-025-57215-9) и [репозиторий авторов](https://github.com/maranasgroup/CatPred) | Контекст базы, OOD-оценка и неопределённость | Локальная копия не идентифицирована с конкретным upstream commit |

Локальные CSV и PDB — отдельный источник фактов о наличии данных;
их хеши и выбранные записи находятся в `evidence/pilot_records.json`.
Сведения о машине получены локальными `lscpu`, `free`, `df` и проверкой
доступных устройств, а не взяты из литературы.

Формулы простой кинетической схемы и оценка чувствительности барьера
получены аналитически в документах. Предлагаемые этапы, ограничения ресурсов
и критерии проверки — рекомендации для этого проекта; они не представлены
как результаты перечисленных статей.
