# Журнал текущего исследования P00947

Новые записи помещаются выше старых; прежние не изменяются. Время локальное MSK.
Записи описывают только действия текущего исследования; более ранние расчёты
имеют собственный журнал [all-events.jsonl](all-events.jsonl).

## 2026-11-02 21:56 | Проверка первого окна с начальным кадром
- В `checked_blocks` исправлен пропуск кадра 0 при энергетической проверке: первое окно теперь 0–4 ps; отдельная температурная проверка, как прежде, использует кадры после начального. Порог и seed не менялись.
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q enzyme-kinetics/tests/test_lowT_quench.py`. Результат 7 passed за 0.42 s; синтетический тест общего генератора останавливается после кадра 4, а не 5. Это не сквозной CLI-тест.
- SHA-256: `lowT_quench.py` 53845ab35a8ecad9c41dfbef19720887f3cc0f6eeb45719ff42a1cc19671577d; `test_lowT_quench.py` dd8481bcf78890f8e6166865ad7f7930201294c7caa28806e76f3d9c46416e94. Новая MD и свип не запускались.
- Блокеры: CLI жёстко требует конкретный лиганд и GLU43 с индексом 601; сквозной тест синтетического пептида отсутствует, summary не имеет полей stop_reason и stop_time_ps. Независимая атомная карта не проверена. Свип не разрешён к запуску до прохождения обоих пользовательских gates.

## 2026-11-02 21:29 | Проверка periodic backbone restraint
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/lowT_quench.py` `a94ef655bb1cfcedb998e80d8eb58939d3309e96040587e7bac5923911b2b8b8`; `enzyme-kinetics/tests/test_lowT_quench.py` `d91d5652f963a6e289e3f32c88929b46c77daba862feb42abf8cec0b9849e600`; прочие тесты проекта. Реальные координаты и XML не использовались.
- Параметры: синтетический OpenMM CPU Context, 4 потока, periodic k*r², k=5 kcal/(mol·Å²), смещение 0.1 нм, тест периодического сдвига 2 нм; новых MD-свипов нет.
- Результат: 83 passed, 1 skipped за 6.11 s; энергия 20.92 кДж/моль и сила −418.4 кДж/(моль·нм) воспроизведены на одноатомной системе. Вывод stdout; summary отсутствует.
- Ограничения: CLI-прогон B, проверка на полной топологии, отдельная запись изменённого System XML и декомпозиция сил пока не выполнены; AllBonds и полный CLI-тест не завершены.
- Статус: partial.

## 2026-11-02 21:27 | Анализ сохранённых исторических кадров GLU43
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/analyze_glu43_history.py --frames enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/quench-298k --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --topology-pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --outdir enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/glu43-history-v1`.
- Входы SHA-256: XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; топология `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; скрипт `420a95738a92e30a0f31d964df907b39400e6617ff472671eb410fa868055872`; полный список SHA-256 всех 21 State/NPY/box и вспомогательного кода — в summary. Граничные State 00 `ef6aa2b28c8d66a39278adf6536952f6c9025e39034d850169ef5cf5445f294a`, 20 `ee43970d9aace27765f2bb29da5fd8908db6b82e1bc4cd37e07547bde488aea0`.
- Параметры: OpenMM 8.6.1, SciPy/NumPy, CPU до 4 потоков; минимум образа для ортогонального box, переразложение сохранённых сил без интегрирования; локальная область 8 Å вокруг backbone C GLU43 A, пары <2 Å без 1–2/1–3. Seed и температура новых расчётов не применимы.
- Выход: `md-pipeline-101/diagnostics/glu43-history-v1/` — geometry.csv, projected-force.csv, environment.csv, local-pbc-pairs.csv, atom-map.csv, summary.json SHA-256 `b5ea82d5a1974e246d703675580e07eb830599396fb8368f8fa68122fe51fa70`.
- Результат: проекция GLU43 A C 475.37 на 0 ps, 2233.81 кДж/(моль·нм) на 6 ps; omega SER42–GLU43 162.41° → −174.28°, GLU43–PRO44 −172.39° → 173.30°. Всего 714 локальных пар-флагов <2 Å за 21 кадр; число не является количеством clashes. Кадры 7–20 явно помечены post_violation.
- Ограничения: не выполнен новый независимый расчёт физических сил, сохранённые силы только повторно спроецированы; атомная идентичность опирается на прежний порядок PDB; локальная выборка не классифицирует химические H-связи; нет сравнения с кристаллом, force-group декомпозиции или вывода о причине. Исторический запуск не является валидным baseline.
- Статус: partial; новая MD и свип не запускались.

## 2026-11-02 21:26 | Тесты анализа исторических кадров GLU43
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/analyze_glu43_history.py` `420a95738a92e30a0f31d964df907b39400e6617ff472671eb410fa868055872`; `enzyme-kinetics/tests/test_analyze_glu43_history.py` `0e20054a49a6d7a1955f29557d9662d00df7531f4f2f60ba6c854e4f53e6aad6` и остальные тесты проекта. Исторические State/NPY в этом тесте не использовались.
- Параметры: Python 3.13, CPU, синтетические геометрии; две новые проверки диэдрала при периодическом сдвиге и отклонения вырожденной геометрии.
- Результат: 82 passed, 1 skipped за 5.90 s. Свип и новая MD не запускались; полноценная проверка анализа на реальных кадрах ещё впереди.
- Статус: partial.

## 2026-11-02 21:21 | Документ требований этапа 2
- Цель: разделить незавершённую техническую работу, пользовательские решения и ресурсы, необходимые для диагностики GLU43 и сравнительного свипа.
- Действия: чтение существующих документов, кода и сохранённых результатов; проверка метаданных окружения, размера сохранённого каталога и временных меток файлов. Новые MD, минимизация, свип, тесты или анализ сил не запускались; существующие научные артефакты не перезаписаны.
- Результат: создан `enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/REQUIREMENTS-TO-UNBLOCK-STAGE2.md`; других новых расчётных артефактов нет. Существующие разрешения на A–E, seed 42 и энергетический порог сохранены, отсутствие завершённой реализации не представлено как отсутствие пользовательского разрешения.
- Ограничения: оценки инженерного и CPU-времени предварительные; протонные состояния остаются модельными гипотезами. Статус этапа 2 не изменён.
- Статус: документ готов; этап 2 остаётся blocked до выполнения технических работ и диагностики.

## 2026-11-02 21:02 | Тесты заготовки restraints и H-связей
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/lowT_quench.py` `7ea15a91464b0956c9f37b00042ae287daa2d1a51e27961544103d412cdd9ffc`; `enzyme-kinetics/tests/test_lowT_quench.py` `e08163952a5fca0cd2daabd1f941b55d5eb405f69891f0a4f7f9ae9a9fc67937`; другие тесты проекта.
- Параметры: Python 3.13 CPU; синтетическая PBC-геометрия H-связи и прежние тесты; реальный комплекс не запускался, summary отсутствует.
- Результат: 80 passed, 1 skipped за 5.80 s. В код добавлена геометрическая оценка N/O контактов и опциональный periodic backbone restraint; химическое типирование акцепторов предварительное.
- Ограничения: полноценный CLI-тест, проверка энергии restraint на реальной системе и AllBonds отсутствуют, B/E не запускались. Тесты не доказывают физической корректности кинетики или топологии.
- Статус: partial.

## 2026-11-02 20:55 | Стартовая геометрия и назначение GLU43
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/inspect_glu43_initial.py --positions-npy enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/positions-nm.npy --topology-pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --crystal-pdb enzyme-kinetics/data/research/p00947/8CHO-assembly.pdb --outdir enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/glu43-initial`.
- Входы SHA-256: NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; topology `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; system `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; 8CHO `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`; исправленный скрипт `bc3b140848152c9876f92d9ca70dc4887f19442792f9fdb8ffa216e1a137cd94`; geometry `db2042bce8535ab1260cd0f01c10bb10f6e981b3d7902528a0b1cffe1a167127`.
- Параметры: OpenMM 8.6.1, без MD; модельная ячейка minimum image, для локальной геометрии исходного кристалла прямые координаты, модели 0/1 сопоставлены цепям A/B.
- Выход: `md-pipeline-101/diagnostics/glu43-initial/summary.json` SHA-256 `bb3d9f6c7f3f5d3944f63d6e8ee3782bacd6167eb34224ae0ea3ccd07871d5f8`.
- Результат: обе цепи GLU без HE1/HE2, заряды CD +0.8054, OE1/OE2 −0.8188 e; это модельное назначение депротонированного глутамата, не экспериментальное доказательство протонирования. Цепь A: CD–OE1/OE2 кристалл 1.2611/1.2532 Å, модель 1.2531/1.2550 Å; угол 122.37° → 121.95°; backbone C–O 1.2405 → 1.2410 Å.
- Ограничения: близость локальной геометрии не исключает локального напряжения; отсутствие выраженного искажения само по себе не установлено статистически, нет оценки среды, силовой декомпозиции или локального PBC-аудита. Свип не выполнялся.
- Статус: partial.

## 2026-11-02 20:54 | Неудачный разбор стартовой геометрии GLU43
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/inspect_glu43_initial.py --positions-npy enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/positions-nm.npy --topology-pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --crystal-pdb enzyme-kinetics/data/research/p00947/8CHO-assembly.pdb --outdir enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/glu43-initial`.
- Входы SHA-256: NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; PDB топология `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; system XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; crystal `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`; скрипт `7ca21676ac23acf14b25c87d3d92fe674466e6e3263dcdf395296ee963fdc9ed`.
- Параметры: Python/OpenMM 8.6.1, CPU, без MD. Результат: exit 1, AttributeError: `getDefaultPeriodicBoxVectors()` вернул список Quantity, а код ожидал один Quantity; каталог и summary не созданы.
- Отрицательный результат: да. Требуется исправить разбор каждого вектора и повторить; предыдущее отрицательное событие сохранить.
- Статус: partial.

## 2026-11-02 20:53 | Проверка остановки в общем цикле NVT
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/lowT_quench.py` `bf7c07b25ae1816eae13b4feb21b72121d0b9eba7051ad0db5997111d4c03086`; `enzyme-kinetics/tests/test_lowT_quench.py` `d796afcffecdaab0dc405252cbc49265855b2ec60758ea952b3b339a6ca720f3` и прочие тесты проекта; физические файлы комплекса не использовались.
- Параметры: Python 3.13, CPU OpenMM, 4 потока для синтетической системы; простой двухчастичный фрагмент и ограниченная трёхатомная вода, искусственный временной энергетический offset −200 кДж/моль на блок, без изменения сил; созданный temporary State XML после каждого блока.
- Результат: 79 passed, 1 skipped за 5.91 s; общий `checked_blocks` остановил цикл до блока 6 на первом проваленном окне 1–5, сохранены synthetic State 0–5; summary.json отсутствует.
- Ограничения: это интеграционный тест общего цикла с реальным Context, но не сквозной запуск основного CLI с метриками и DCD; полноценная валидация CLI остаётся открытой. Зависит от искусственного потенциала, не модели KSI.
- Статус: partial.

## 2026-11-02 08:35 | Исправление и тест скользящей проверки
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: исправленный `enzyme-kinetics/scripts/lowT_quench.py` `e118683ce725caf9d668525612691fc187afe891f07c31a9490c8c9e7fe64d11`; `enzyme-kinetics/tests/test_lowT_quench.py` `a473181e512beddfc9a34c4be08adec941d636c10e5e7d40079ac5ee07600d70`; прочие тесты. Числа из сохранённого metrics.csv `d4f79efd6a7b1d901805343af285f6f0936620e59f7eaefe2826d9a9b8e13f2f` воспроизведены в тесте, файл самим тестом не читается.
- Параметры: Python 3.13, CPU; проверка функций и регрессий без новой динамики; выход stdout, summary.json отсутствует.
- Результат: 78 passed, 1 skipped, 5.60 s. Условие energy/temperature теперь проверяется после каждого сохранённого блока при наличии пяти кадров; unit-тест подтверждает первое нарушенное окно, но не является сквозным тестом остановки MD.
- Отрицательный результат?: исходный прогон уже нарушил требование ранней остановки и не исправляется изменением кода; повтор не проводился.
- Статус: partial, этап 2 blocked.

## 2026-11-02 08:35 | Исправление интерпретации проверки энергии NVT
- Цель: проверить условие остановки на каждом окне пяти сохранённых кадров без повторного MD.
- Команда: чтение `md-pipeline-101/quench-298k/metrics.csv` через `csv.DictReader`, поиск строго монотонного падения >500 кДж/моль среди скользящих окон пяти кадров 1–20 ps.
- Входы: сохранённый metrics.csv из запуска 08:34, без новых сил или координат; остальные хеши входов вычисления в предыдущей записи.
- Результат: первое нарушающее окно 2–6 ps, ΔE = −971.73 кДж/моль; последующие окна 3–7, 15–19 и 16–20 ps также не прошли. Скрипт проверил критерий только после последнего блока и потому ошибочно продолжил 14 ps после требуемой остановки.
- Отрицательный результат?: да — нарушено требование немедленно остановиться на первом нарушении; данные 7–20 ps сохранены как история этого ошибочно продолженного запуска, не использовать как результат прошедшего протокол этапа. State 6 ps сохранён, но никакой автоматический restart не разрешается.
- Влияние на следующие шаги: исправить и протестировать покадровую скользящую проверку; не выполнять повтор или эквилибрацию без согласования. Температура и физический смысл потенциальной энергии требуют отдельной оценки; не объявлять готовность к MD.
- Статус: blocked.

## 2026-11-02 08:34 | Short NVT MD test при 298.15 K: отрицательный результат
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/lowT_quench.py --positions-npy enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/positions-nm.npy --state-xml enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/state.xml --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --box-vectors enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/box-vectors-nm.npy --restraints-json enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/restraints.json --topology-pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --temperature 298.15 --friction 10.0 --timestep-fs 1.0 --duration-ps 20 --save-interval-ps 1 --seed 101 --outdir enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/quench-298k`.
- Входы SHA-256: NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; State `ba90694d649da362333255dad135fd3469e0797d480814d0ed3b7b2188242cf1`; box `224b194b7e1ca1c225f5167fcbce1faec9e7b31bd79ac056d6acc2db4baa4342`; system XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; restraints `4392227bc0e43cb421c99b27827a7908d45e5a2d36fb12b46a69b1c89258c919`; PDB топология `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; скрипт `fb3a136d35fdc60bbeeb088d8c1b1592898073fe9f1bc9fecd1443bc8ca593ad`; остальные хеши в summary.
- Параметры: OpenMM 8.6.1 CPU 4 потока, NVT LangevinMiddle 298.15 K, 10 ps⁻¹, 1 fs, 20 ps, кадры/State каждую 1 ps, seed 101; 41 515 constraints из XML, HBonds при подготовке, tolerance интегратора 1e-5 относительный, мягких pose restraints нет; начальные скорости сгенерированы заново; restart from coordinates, not bit-exact continuation.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/quench-298k/`, summary SHA-256 `d0e95a757799ca6cc608bad29e4c9d87c9a60e5b83897b386c581b12af8f5252`; run-config SHA-256 `ab7bcb0a3c7c7abbb9e778b12e180d8cdbb110436a8a5a6f42769fa92cdad971`; 21 сохранённый State включая старт, 20 кадров DCD.
- Результат: 20 ps рассчитаны, но код выхода 1 и status failed: последние пять Epot монотонно −553587.59 → −554985.83 кДж/моль, изменение −1398.24 (>500). Средняя температура последних пяти 298.65 K; финальная 299.97 K; финальные RMSD лиганда 0.07067 нм, backbone 0.07546 нм, контакты 15; финальная ошибка constraints 7.48e-8 нм. GLU43 A backbone C Fproj 475.37 → 1372.24 кДж/(моль·нм).
- Отрицательный результат?: да — заданный критерий энергии не пройден; не считать систему готовой к эквилибрации и не продолжать без согласования. Мгновенные силы в 298 K флуктуируют и не сравниваются напрямую с минимизацией как доказательство ухудшения.
- Ограничения: отсутствует полный интеграционный тест NVT-скрипта на синтетическом пептиде с водой; геометрические H-связи кандидаты, атомные индексы топологии предполагают прежний порядок; короткая MD не доказывает равновесия или WT–ES.
- Статус: blocked; нагрев, production MD, enhanced sampling и QM/MM не запускались.

## 2026-11-02 08:15 | Тесты проверки short NVT
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=4 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/lowT_quench.py` `fb3a136d35fdc60bbeeb088d8c1b1592898073fe9f1bc9fecd1443bc8ca593ad`; `enzyme-kinetics/tests/test_lowT_quench.py` `6420bd291ba8b296580a5edf84b5273360458d77a3deb06c64ab14322cc5264e`; остальные тесты проекта. Реальная молекулярная система не использовалась.
- Параметры: Python 3.13, CPU; синтетические данные и сериализация одноатомного State; ограниченная проверка защитных функций, не полная интеграция MD на модельном пептиде.
- Результат: 77 passed, 1 skipped за 5.75 s; stdout, summary.json отсутствует.
- Ограничения: нет теста полного NVT-контура на пептиде и воде; результаты тестов не гарантируют стабильности комплекса и не являются научной валидацией.
- Статус: partial.

## 2026-11-02 08:13 | Тест геометрии GLU43 для планируемого NVT
- Цель: отделить геометрию backbone C с высокой силой от карбоксилата боковой цепи и проверить periodic minimum image.
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=1 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы: скрипт `enzyme-kinetics/scripts/glu43_geometry.py` SHA-256 `db2042bce8535ab1260cd0f01c10bb10f6e981b3d7902528a0b1cffe1a167127`; тест `enzyme-kinetics/tests/test_glu43_geometry.py` SHA-256 `92f536e897f509eb3e489ce4e4fbee13e5c3bbe1f8c9d55b91262fafb9349ff8`; прочие тесты проекта. NPY, State, box, system.xml и restraints в тесте не использовались.
- Параметры: Python 3.13, синтетическая геометрия, CPU; симуляция не выполнялась. Вывод: stdout pytest; summary.json отсутствует.
- Результат: 74 passed, 1 skipped, 5.53 s. Тестируется различение backbone C–O и боковой цепи CD–OE1/OE2 и периодический перенос.
- Ограничения: atom C с высокой силой — backbone carbonyl C, а не CD карбоксилата; у backbone C нет пары OE1/OE2. Проверка не является запуском NVT или диагностикой реальной временной эволюции, интеграция в NVT-скрипт ещё не выполнена.
- Статус: partial; `lowT_quench.py` и MD на 298.15 K не выполнены.

## 2026-11-02 07:54 | PBC-аудит контактов модели 101
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/pbc_clash_audit.py --input-npy enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/positions-nm.npy --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --topology-pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --state-xml enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/state.xml --box-npy enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/box-vectors-nm.npy --restraints-json enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force/restraints.json --outdir enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/pbc-clash`.
- Входы SHA-256: NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; system XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; topology PDB `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; state XML `ba90694d649da362333255dad135fd3469e0797d480814d0ed3b7b2188242cf1`; box `224b194b7e1ca1c225f5167fcbce1faec9e7b31bd79ac056d6acc2db4baa4342`; restraints `4392227bc0e43cb421c99b27827a7908d45e5a2d36fb12b46a69b1c89258c919`; script `7f15cef8c4d9d291e8e102192573fdffca1ea34fb12bb97d9409d881d692c594`.
- Параметры: OpenMM 8.6.1, CPU, ортогональная ячейка 7.6980177893 нм, minimum image, интервалы <1.5 и <2 Å, категории и исключения указаны в summary; только чтение State, динамики нет.
- Выходы: `md-pipeline-101/diagnostics/pbc-clash/`, summary SHA-256 `d63343052aafa2ad357ce41f337995adaeadd8704cfaeb7e617dbcc095c182d8`, state-used SHA-256 `c6453d6b54f2f0fd4d36397a4a0a45dbd4599e96c5e01bd1db97dc333e4c3e49`.
- Результат: <1.5 Å: 0 во всех выбранных категориях; <2 Å: вода–вода 23465, белок–белок 12, лиганд–белок 0, лиганд–вода 0. Все перечисленные пары содержат H; 386 молекул пересекают границу wrapped box по связанным атомам.
- Ограничения: большинство коротких O–H пар воды может соответствовать водородным связям, порог 2 Å сам по себе не доказывает физический clash. Только указанные межмолекулярные категории, PDB только топология при предположении сохранённого порядка атомов; не глобальный всепарный аудит.
- Статус: done для указанной ограниченной проверки; quench ещё не запущен.

## 2026-11-02 07:54 | Тесты PBC-аудита
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=1 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/pbc_clash_audit.py` `7f15cef8c4d9d291e8e102192573fdffca1ea34fb12bb97d9409d881d692c594`; `enzyme-kinetics/tests/test_pbc_clash_audit.py` `f91e6f6a204bdcf773f54e6be45c8e9b42c644e7c17c17e26f8c5801baf9dee4`; прочие тесты проекта. Молекулярные NPY, State и system XML не использовались.
- Параметры: Python 3.13, OMP_NUM_THREADS=1, синтетические данные; OpenMM-расчёт не проводился.
- Результат: 72 passed, 1 skipped, 5.51 s. Проверены minimum image, сдвиг на ячейку и логика исключений. stdout; summary.json отсутствует.
- Ограничения: тесты не являются полным доказательством корректности атомного маппинга и химической интерпретации порогов.
- Статус: done для тестов.

## 2026-11-01 22:07 | Тесты диагностики MD pipeline
- Цель: проверить группировку сил и регрессии проекта.
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=1 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы: `enzyme-kinetics/tests/test_diagnose_high_force_atoms.py` SHA-256 `e391b4d160c03255b40c60fc55f533c7e9b22a6b51db41020e15e5f3542820e6`; скрипт `717a3bf2c8b13f76a6184a89ab6776ed92d830d74552b32e739a51d7e84c65c2`; остальные тесты проекта.
- Параметры: Python 3.13, CPU, OMP_NUM_THREADS=1; синтетическое распределение по категориям. Молекулярные входы не использованы.
- Результат: 70 passed, 1 skipped, 5.52 s; stdout pytest, summary.json не создавался.
- Ограничения: тест не проверяет корректность всего маппинга или физическую готовность системы; PBC clash и quench не выполнены.
- Статус: partial.

## 2026-11-01 22:06 | MD pipeline: распределение проецированных сил
- Цель: локализовать атомы с высокой проецированной силой без запуска динамики.
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/diagnose_high_force_atoms.py --input-npy enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-continuation/stage-2-positions-nm.npy --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --topology-pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --restraints-json enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-convergence-v2/restraints-reconstructed.json --outdir enzyme-kinetics/outputs/active-site-v3/md-pipeline-101/diagnostics/high-force`.
- Входы SHA-256: NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; PDB только для топологии `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; restraints `4392227bc0e43cb421c99b27827a7908d45e5a2d36fb12b46a69b1c89258c919`; скрипт `717a3bf2c8b13f76a6184a89ab6776ed92d830d74552b32e739a51d7e84c65c2`; projection `b7814919d383ba62b60330b5db00d03f8580ed369fc4667146950899d3000c85`.
- Параметры: OpenMM 8.6.1, CPU 4 потока, без оптимизации и динамики; порог высокой силы 50 кДж/(моль·нм), 41 515 constraints, координаты NPY, cubic box из system XML.
- Выходы: `md-pipeline-101/diagnostics/high-force/`: per-atom.csv, top100.csv, state.xml, box-vectors-nm.npy, positions-nm.npy, restraints.json, summary.json SHA-256 `1d7c7c6d272b77e1267223ad21b1de14c5b5f509665b80998c83b032f68bec5a`.
- Результат: 42 атома >50; белок 35 (83,3%), вода 7 (16,7%), лиганд/ионы 0. Максимумы белок 475,37, лиганд 24,88, вода 108,98, ионы 15,64 кДж/(моль·нм); первый атом GLU43 A C. Остаток проекции 1,77e-8.
- Отрицательный результат?: PBC clash аудит ещё не выполнен, высокий максимум белка требует дополнительной диагностики; пригодность к MD и происхождение силы неизвестны. Маппинг топологии опирается на сохранённый порядок PDB.
- Влияние на следующие шаги: не запускать quench до проверки PBC контактов и реализации проверок; этапы 3–6 не запускались.
- Статус: partial.

## 2026-11-01 21:47 | Проекция сил полной системы на фиксированных координатах
- Цель: проверить численную ортогональность проекции и точность ограничений без минимизации.
- Команда: `OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 "$HOME/.local/share/enz-kinetics-mm/bin/python" enzyme-kinetics/scripts/check_projection_full_system.py --input-npy enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-continuation/stage-2-positions-nm.npy --system-xml enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --output enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-convergence-v2/projection-check/summary.json`.
- Входы SHA-256: NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; system XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; скрипт `b77ad27cf8ca031a2ace7dee0acea4655c4e58dffb0775d233c27c78c1644563`; проекция `b7814919d383ba62b60330b5db00d03f8580ed369fc4667146950899d3000c85`. Реконструированный restraints JSON `4392227bc0e43cb421c99b27827a7908d45e5a2d36fb12b46a69b1c89258c919` не использован этим запуском.
- Параметры: OpenMM 8.6.1 CPU, 4 потока, фиксированные NPY координаты, 41 515 ограничений, допуск интегратора 1e-5 относительный; блоки, seed и max_iterations не применимы, минимизация не выполнялась.
- Выход: `enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-convergence-v2/projection-check/summary.json`, SHA-256 `8e8aa58c26aa46e654021369dea7519f1621d2ebd4d958d085ca4d612d62fcf6`.
- Результат: max |Fproj| 475.3671 кДж/(моль·нм), RMS компонентов 3.88323; max |J Fproj| 1.77e-8 кДж/(моль·нм); максимальная ошибка длины 4.67e-9 нм, RMS 2.67e-9 нм, max относительная 3.085e-8, сверх относительного допуска 0. Проекция меняет силы, не координаты, поэтому геометрические ошибки до/после тождественны и их неизменность не является независимой проверкой проекции.
- Отрицательный результат?: да — критерий max Fproj < 10 не выполнен в этой точке, и трёх блоков нет. Обычный L-BFGS-B с подменой истинного градиента его проекцией не сохраняет ограничения вдоль конечных шагов и нарушает предпосылки линейного поиска; success SciPy не доказывает сходимость ограниченной задачи.
- Влияние на следующие шаги: до расчёта необходим корректный ограниченный оптимизатор с допустимыми пробными шагами, проверкой геометрии и достоверным stop_reason; PBC/водный аудит и MD не выполнялись.
- Статус: partial.

## 2026-11-01 21:42 | Проверка проекции и критериев
- Цель: проверить математическую основу проекции и строгий критерий трёх блоков, а также регрессии проекта.
- Команда: `PYTHONPATH=.venv/lib/python3.13/site-packages:enzyme-kinetics/src OMP_NUM_THREADS=1 /home/vladmyahlov/.local/share/uv/python/cpython-3.13-linux-x86_64-gnu/bin/python3.13 -m pytest -q`.
- Входы SHA-256: `enzyme-kinetics/scripts/constraint_projection.py` `b7814919d383ba62b60330b5db00d03f8580ed369fc4667146950899d3000c85`; `enzyme-kinetics/tests/test_constraint_projection.py` `67383347bdb843b443bf5c5be56662104fab448ab22e07a8737a7d97a744ee83`; остальные тесты проекта. Молекулярные входы NPY/XML/restraints не использовались; summary.json отсутствует.
- Параметры: Python 3.13, OMP_NUM_THREADS=1; синтетическая задача без OpenMM. LSMR, проверка остатка J·Fproj; блок 100, три последовательных блока, заданные пороги без изменений.
- Результат: 69 passed, 1 skipped за 5.50 s. Проверены аналитическая пара, зависимые ограничения, вырожденная геометрия и критерий плато.
- Отрицательный результат?: да — проекция пока не валидирована на полной системе, минимизатор и PBC/водный аудит не реализованы, перезапуск не выполнялся. OpenMM MinimizationReporter сообщает итерации штрафной оптимизации, которые могут повторяться или отменяться, и не возвращает причину остановки.
- Влияние на следующие шаги: необходима проверенная инструментированная оптимизация с достоверным stop_reason и сохранением блоков; MD и QM/MM не запускались.
- Статус: partial.

## 2026-11-01 21:40 | Предрегистрация и реконструкция входов v2
- Цель: зафиксировать критерии и проверить источники restart после разрешения пользователя использовать кубическую simulation box из system.xml, несмотря на несовпадение с кристаллографической CRYST1.
- Команда: только чтение и разбор `system.xml` через Python `xml.etree.ElementTree`, просмотр `CRYST1` и скриптов подготовки/аудита; расчёт OpenMM и тесты не запускались.
- Входы SHA-256: system XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; конечный NPY `011449cf2769150fbd5982095da247ff6b288bbcd214a95e1817a831032507da`; 8CHO `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`; audit script `f42888dc8e6708d7928deff9566a32627458c973087577b002d17da405198d4c`; preparation script `cb17d87a8c99df973fb0656b4301e0ec1206bf813ff359388d05749981f92e0e`.
- Параметры: только предрегистрация; план CPU 4 потока, блок 100, три последовательных блока с |ΔE| < 0,1 кДж/моль, max Fproj < 10 кДж/(моль·нм), max Δr < 0,01 Å и stop_reason=converged. seed, max_iterations и tolerance нового расчёта ещё не заданы. Мягкие pose restraints не обнаружены; 41 515 жёстких distance constraints сохранены в system XML.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-convergence-v2/` — `convergence-config.json`, `restraints-reconstructed.json`, предварительный `REPORT.md`; summary.json отсутствует.
- Результат: finite cubic box 7,698017789305977 нм, баростат в XML отсутствует; CRYST1 6,168 × 6,168 × 14,360 нм, 90°/90°/120°. Согласовано пользователем, что применяется simulation box, а не CRYST1. Классификация: instrumented restart from coordinates, not bit-exact continuation; State XML предыдущего расчёта не сохранён.
- Отрицательный результат?: да — новый минимизатор и аудит не выполнены; сходимость, стоп-причина, проекция и PBC-контакты неизвестны. Отличие box от кристалла не является само по себе доказательством ошибки.
- Влияние на следующие шаги: реализовать и проверить инструментированный расчёт, выполнить PBC/водный аудит; MD и QM/MM не запускались.
- Статус: partial.

## 2026-11-01 20:38 | Тесты после точного restart
- Цель: проверить регрессии после добавления restart из NPY.
- Скрипт/команда: `python -m pytest -q` при `OMP_NUM_THREADS=1`.
- Входы: скрипт `enzyme-kinetics/scripts/audit_ksi_minimization.py` SHA-256 `f42888dc8e6708d7928deff9566a32627458c973087577b002d17da405198d4c` и тесты проекта.
- Параметры: Python 3.13, проектные site-packages, CPU.
- Выходы: stdout pytest.
- Результат: 66 passed, 1 skipped за 5.84 s; restart физической системы отдельно проверен успешным запуском.
- Отрицательный результат?: нет; визуализационный тест пропущен в core из-за отсутствия imageio и ранее прошёл в tools-окружении.
- Диагностика: unit-тесты не проверяют сходимость или физическую корректность системы.
- Влияние на следующие шаги: функциональных регрессий не обнаружено, критерии научной приёмки остаются открытыми.
- Статус: done.

## 2026-11-01 20:37 | Продолжение с координатами полной точности
- Цель: оценить, прекращается ли снижение энергии после двух дополнительных блоков минимизации.
- Скрипт/команда: `audit_ksi_minimization.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --assignment enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json --positions-npy enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-crystalwater/stage-2-positions-nm.npy --output enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-continuation --iterations 300 --tolerance 1.0 --threads 4`.
- Входы: PDB SHA-256 `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; assignment `794d5149e71874cee9093a1af456559eb7c8c33abdea4a7b6f150a28ae200db8`; NPY `c55878a4a74f9e9fae3eb8ce47970c429c6acd97d680e3935602d484905a22a6`; новая версия скрипта `f42888dc8e6708d7928deff9566a32627458c973087577b002d17da405198d4c`.
- Параметры: OpenMM 8.6.1, CPU 4 потока, 2×300 итераций, tolerance 1 кДж/(моль·нм), атомный порядок и box прежние.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-continuation/`.
- Результат: энергия −678946.91 → −684898.73 → −688074.49 кДж/моль; второй блок ещё −3175.76. OD1 C4–O 3.558 → 3.622 Å, Hβ–O 2.500 → 2.560 Å, угол 161.50 → 162.71°.
- Отрицательный результат?: да — уменьшение энергии сохраняется, заявленная сходимость не установлена; число выполненных итераций и причина остановки недоступны в данном выводе.
- Диагностика: старт по NPY воспроизводит предыдущую энергию с разницей 0.0047 кДж/моль из-за численного вычисления, без нового PDB-округления; сырые силы содержат constraint-компоненты.
- Влияние на следующие шаги: модель всё ещё нельзя объявлять сошедшейся или равновесной; требуется приборный критерий сходимости и проверка сил/ограничений, а не бесконечное увеличение итераций.
- Статус: partial.

## 2026-11-01 20:32 | Тесты аудита минимизации
- Цель: проверить новый код и регрессии проекта.
- Скрипт/команда: `python -m pytest -q` при `OMP_NUM_THREADS=1`.
- Входы: `enzyme-kinetics/scripts/audit_ksi_minimization.py` SHA-256 `2518b546dd5cdb908a416d440e6ef104db24f3cc1f31d1696bb3d6f09151423a`; `enzyme-kinetics/tests/test_audit_ksi_minimization.py` SHA-256 `9b26084d69862ea780bbf7048819ba71d3a4bcf528e20b5ab4427e70ce12a10b` и прочие тесты проекта.
- Параметры: Python 3.13, проектные site-packages, CPU.
- Выходы: stdout pytest.
- Результат: 66 passed, 1 skipped за 5.63 s.
- Отрицательный результат?: нет; необязательный imageio-тест пропущен в этом окружении, отдельно выполнялся в tools-окружении ранее.
- Диагностика: добавлены тесты корректности вычисления геометрии, сил и отклонения NaN; тесты не доказывают сходимость физической системы.
- Влияние на следующие шаги: программные регрессии не обнаружены.
- Статус: done.

## 2026-11-01 20:31 | Аудит и продолжение минимизации KSI
- Цель: проверить исходную физическую модель 101 с кристаллографическими водами и приблизиться к сходимости минимизации.
- Скрипт/команда: `audit_ksi_minimization.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --assignment enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json --output enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-crystalwater --iterations 300 --tolerance 1.0 --threads 4`.
- Входы: PDB SHA-256 `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`; XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`; assignment `794d5149e71874cee9093a1af456559eb7c8c33abdea4a7b6f150a28ae200db8`; скрипт SHA-256 `2518b546dd5cdb908a416d440e6ef104db24f3cc1f31d1696bb3d6f09151423a`.
- Параметры: OpenMM 8.6.1 CPU, 4 потока, 2 последовательных блока до 300 итераций, tolerance 1 кДж/(моль·нм), 43 434 атома, 41 515 ограничений.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/minimization-audit-101-crystalwater/summary.json`, PDB и точные массивы координат для трёх стадий.
- Результат: энергия −616065.13 → −668074.66 → −678946.92 кДж/моль; второй блок снизил её на 10872.25. Суммарный заряд ≈0 e, заряд лиганда ≈0 e; нечисловых параметров nonbonded и отрицательных σ/ε нет. C4–OD1 3.520 → 3.558 Å, Hβ–OD1 2.461 → 2.500 Å, угол 161.69 → 161.50°. 48 исходных вод прослежены, сдвиг O относительно начала аудита: медиана 0.680 Å, максимум 2.119 Å. Лиганд–белок тяжёлые атомы ближе 2 Å не обнаружены.
- Отрицательный результат?: да — строгая сходимость не доказана: энергия за второй блок существенно изменилась; выход минимизатора не содержит причины остановки. Сырые силы включают силы геометрических ограничений и не являются проекционным критерием остановки.
- Диагностика: перезапуск из округлённого PDB; проверка столкновений только тяжёлые атомы лиганда против белка без PBC, не полный clash-аудит. Энергия включает всю систему/воду и не является энергией связывания; смещение вод относительно начала аудита, не исходного кристалла.
- Влияние на следующие шаги: не запускать как доказанный равновесный ES; необходимы сходимость по ограниченной оптимизации, расширенный clash-анализ и контролируемая подготовка альтернативных микросостояний.
- Статус: partial — вычисления завершены, критерий сходимости не выполнен.

## 2026-11-01 13:49 | Обновлённая фиксация tools-окружения
- Цель: сохранить обновлённый состав tools после установки pytest.
- Скрипт/команда: `uv pip freeze --python "$HOME/.local/share/enzyme-kinetics-tools/bin/python" > docs/experiments/active-site-v3/tools-pip-freeze-v2.txt`.
- Входы: tools-окружение Python 3.12.
- Параметры: uv pip freeze, без изменения установленных пакетов.
- Выходы: `docs/experiments/active-site-v3/tools-pip-freeze-v2.txt`.
- Результат: новый снимок содержит pytest 9.2.0; первоначальный снимок не изменён.
- Отрицательный результат?: нет.
- Диагностика: версии не гарантируют неизменность источников пакетов.
- Влияние на следующие шаги: можно воспроизводить отдельный тест визуализации.
- Статус: done.

## 2026-11-01 13:48 | Тест рендера в tools-окружении
- Цель: выполнить пропущенный из core тест визуализации.
- Скрипт/команда: `uv pip install --python "$HOME/.local/share/enzyme-kinetics-tools/bin/python" pytest`; затем `python -m pytest -q enzyme-kinetics/tests/test_render_ksi_preview.py` в том же окружении.
- Входы: `enzyme-kinetics/tests/test_render_ksi_preview.py`, tools-окружение Python 3.12; PyPI pytest.
- Параметры: pytest 9.2.0, imageio установлен.
- Выходы: вывод pytest; обновлено tools-окружение.
- Результат: 1 passed за 0.60 s.
- Отрицательный результат?: нет; этот unit-тест не проверяет визуальную точность изображения.
- Диагностика: ранее этот тест пропущен только в проектном окружении без imageio; само изображение создано отдельным функциональным запуском.
- Влияние на следующие шаги: все добавленные тесты выполнены, но в разных изолированных окружениях; список tools-версий требует обновления.
- Статус: done.

## 2026-11-01 13:46 | Исправление хеша из записи 13:45
- Цель: исправить ошибку транскрипции SHA-256.
- Скрипт/команда: сверка вывода `sha256sum enzyme-kinetics/scripts/check_qmmm_embedding.py`.
- Входы: `enzyme-kinetics/scripts/check_qmmm_embedding.py`.
- Параметры: SHA-256.
- Выходы: эта запись; предыдущая запись сохранена без изменений.
- Результат: правильный SHA-256 `21a6eac07fbd57ba21bf052b551d16d0d014b9ec422c011fe1e8c9813e31e787`.
- Отрицательный результат?: да — в записи 13:45 для этого файла ошибочный хеш.
- Диагностика: остальные перечисленные хеши перепроверять при воспроизведении.
- Влияние на следующие шаги: использовать исправленное значение.
- Статус: done.

## 2026-11-01 13:45 | Хеши диагностических скриптов
- Цель: зафиксировать версии кода, связанного с текущими результатами.
- Скрипт/команда: `sha256sum` для файлов в `enzyme-kinetics/scripts/`.
- Входы: `prepare_ksi_complex.py` `cb17d87a8c99df973fb0656b4301e0ec1206bf813ff359388d05749981f92e0e`; `run-rusticl-openmm.sh` `c5b679b09ad8dda344859e332bdd2d0b0144071b2fa2ef6c564d111301f526e5`; `validate_opencl_gpu.py` `ace5b4113913949c7997eaba0d731bb4b51ace7a0556fa91cb9e6b625b1ae112`; `verify_ksi_atom_mapping.py` `ccd6e76fe801c0dd46d5871fb2878c0ed83ee0cc4d813e357adf2b0cc1bfd17a`; `check_qmmm_embedding.py` `21a6eac07fbd57ba45a08a7e4011194e1fd6b074090bdb35794c5df3c7fb5f22`; `render_ksi_preview.py` `73166e7c8d7153c0446a8be6a7ca7051e0e1e69f644c7527e77d19cdab75ffb0`.
- Параметры: SHA-256 текущих файлов.
- Выходы: эта запись.
- Результат: хеши записаны.
- Отрицательный результат?: нет.
- Диагностика: для отдельных предыдущих запусков коды могли отличаться; ссылка на хеш не заменяет проверку версий на момент запуска.
- Влияние на следующие шаги: улучшение воспроизводимости.
- Статус: done.

## 2026-11-01 13:44 | OpenCL: возможности Radeon
- Цель: объяснить недоступность mixed precision на текущей платформе.
- Скрипт/команда: локальный `clinfo` через bwrap с Rusticl/radeonsi; фильтр сведений о GPU и double precision.
- Входы: пользовательский Mesa OpenCL и GPU AMD `[1002:699f]`; состояние устройств динамическое.
- Параметры: Rusticl 25.2.8, OpenCL 3.0, Polaris12.
- Выходы: stdout.
- Результат: RX 550 доступен как GPU, глобальная память 8 GiB, double-precision floating point обозначена `(n/a)`.
- Отрицательный результат?: да — нужная двойная точность не доступна через текущий OpenCL runtime; это согласуется с отказом mixed, но не доказывает причину segfault полного single расчёта.
- Диагностика: PCI указывает общее семейство, уточнённое имя устройства дано драйвером; оценка памяти не является измерением свободного объёма.
- Влияние на следующие шаги: single-only GPU остаётся неподдерживаемым OpenMM для production, CPU — корректный fallback.
- Статус: partial.

## 2026-11-01 13:43 | Тесты после добавления кристаллографической воды
- Цель: проверить регрессии проекта после новых расчётов и изменения подготовки комплекса.
- Скрипт/команда: `python -m pytest -q` с OMP_NUM_THREADS=1.
- Входы: `enzyme-kinetics/tests/` и код проекта, включая изменённый `enzyme-kinetics/scripts/prepare_ksi_complex.py`; отдельный хеш изменённого скрипта фиксируется в последующей записи.
- Параметры: Python 3.13 и существующий site-packages, CPU.
- Выходы: вывод pytest.
- Результат: 64 passed, 1 skipped за 5.50 s.
- Отрицательный результат?: нет; тест рендера пропущен из-за отсутствия imageio в этом окружении; тесты сами по себе не проверяют корректность QM/MM и стабильность GPU.
- Диагностика: функциональный рендер и минимизация выполнялись в отдельных окружениях ранее.
- Влияние на следующие шаги: кодовая регрессия не обнаружена; требуется научная валидация.
- Статус: done.

## 2026-11-01 13:40 | Смещение кристаллографических кислородов
- Цель: проверить, остались ли 48 заданных O воды в построенной структуре после минимизации.
- Скрипт/команда: Bio.PDB/numpy; последовательное сравнение вод обоих моделей `8CHO-assembly.pdb` с O остатков HOH цепи D в минимизированном PDB.
- Входы: сборка 8CHO SHA-256 `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`, результат PDB `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`.
- Параметры: сравнение в общей системе координат, соответствие по сохранённому порядку без выравнивания.
- Выходы: stdout.
- Результат: 48/48 атомов O; медианное смещение 0.360 Å, максимальное 0.852 Å, выше 1 Å нет.
- Отрицательный результат?: нет; сохранение положения при 100 итерациях не доказывает устойчивость, каталитическую роль или достоверность водородов воды.
- Диагностика: сохранение порядка соответствует способу построения; проверка не включает структурные воды из других кристаллических моделей.
- Влияние на следующие шаги: можно анализировать сохранённую водную сеть как модельный старт, но требуются равновесные выборки.
- Статус: done для проверенных позиций, partial для биофизической валидации.

## 2026-11-01 13:39 | Геометрия комплекса с водой 8CHO
- Цель: проверить β-H после минимизации с 48 кристаллографическими O воды.
- Скрипт/команда: `analyze_ksi_complex.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/system.xml --beta-assignment enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/geometry-interaction.json`.
- Входы: PDB SHA-256 `77e50f95e2200ab9cce88570eadc61863b8e01ff3b515b6642358148b16bbb4e`, XML `26ceb556460b25b04b0595e138a95f77a758d046d139db85db218e97148a8201`, assignment `794d5149e71874cee9093a1af456559eb7c8c33abdea4a7b6f150a28ae200db8`.
- Параметры: прямой непериодический Coulomb+LJ белок–лиганд; геометрия OD1 и OD2.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/geometry-interaction.json`.
- Результат: OD1: C4–O 3.520 Å, Hβ–O 2.461 Å, угол 161.69°; прямой белок–лиганд −142.17 кДж/моль.
- Отрицательный результат?: нет; энергия не является ΔG связывания, вода исключена из прямой суммы и 100 итераций не доказывают сходимость.
- Диагностика: сравнение с более ранним комплексом даёт похожую геометрию, но системы и молекулы воды различны, статистика отсутствует.
- Влияние на следующие шаги: измерить смещения исходных вод и провести независимые проверки релаксации и протонирования.
- Статус: partial.

## 2026-11-01 13:37 | Комплекс 101 с исходными кислородами вод
- Цель: сохранить 48 кислородов кристаллографической воды 8CHO в минимизации.
- Скрипт/команда: `prepare_ksi_complex.py --protein enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb --ligand-sdf enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/ligand.sdf --crystal-assembly enzyme-kinetics/data/research/p00947/8CHO-assembly.pdb --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater --max-iterations 100 --padding-nm 0.8`.
- Входы: WT SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`, SDF `752c8bcbf25d6fd60bbde6696d0183fc9de647e89f2970e4304fd97844725b19`, сборка 8CHO `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`.
- Параметры: Amber14/GAFF2.11/AM1-BCC, TIP3P, pH 7, ASH99 обеих цепей, 100 итераций, CPU 2; водороды кристаллических вод заданы в произвольной начальной ориентации.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101-crystalwater/`.
- Результат: 48 позиций O воды включены как исходные, 43 434 атома; полная энергия системы −314920 → −616066 кДж/моль.
- Отрицательный результат?: нет, но воды не доказаны каталитическими, их водороды и состояния ориентированы модельно; общая энергия не является энергией связывания.
- Диагностика: скрипт проверяет совпадение 250 CA и отсутствие перекрытия кристаллического O с лигандом менее 2.4 Å; сохранность вод после минимизации и сходимость требуют количественной проверки.
- Влияние на следующие шаги: сравнить координаты вод до/после и геометрию лиганда; не считать сеть равновесной.
- Статус: partial.

## 2026-11-01 13:35 | Перекрытия кристаллографических вод до релаксации
- Цель: проверить, можно ли рассматривать сохранение вод 8CHO в исходных кандидатах.
- Скрипт/команда: Bio.PDB + numpy, расстояние 48 атомов O воды `8CHO-assembly.pdb` до тяжёлых атомов `candidate-{101,202}.xyz`.
- Входы: 8CHO сборка SHA-256 `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`, XYZ 101 `d454619a37b907ac61150f95a38f753b82c9b494b8cd9f2e9d67e35fc01cd862`, XYZ 202 `8c91372c9b301d1dc4b640553a00106f3a87632590f3a4d44c073bfa80dfb482`.
- Параметры: только геометрия исходных фиксированных вод, тяжёлые атомы, пороги 2.0/2.4 Å.
- Выходы: stdout.
- Результат: из 48 вод ни одна не ближе 2.4 Å к тяжёлым атомам стероида; ближайшие 3.046 Å (101) и 2.698 Å (202).
- Отрицательный результат?: нет; сохранность и функциональное значение вод после солватации и минимизации не доказаны.
- Диагностика: модель без водородов воды; критерий перекрытий не доказывает оптимальную водную сеть.
- Влияние на следующие шаги: сохранение исходных вод возможно рассмотреть как отдельную воспроизводимую модель, не автоматически объявляя воды каталитическими.
- Статус: partial.

## 2026-11-01 13:34 | Геометрия позы 202 после минимизации
- Цель: сопоставить положение модельного β-H со второй минимизированной структурой.
- Скрипт/команда: `analyze_ksi_complex.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-202/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-202/system.xml --beta-assignment enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-202/geometry-interaction.json`.
- Входы: PDB SHA-256 `8ae66af312de811c4e5d57c1be026251b08456d5cd6039180e53558dc00b0935`; XML `8842dca400bb459770d8ec7db32c6bde1a43abfb1c4145587056cb9601afb047`; assignment `794d5149e71874cee9093a1af456559eb7c8c33abdea4a7b6f150a28ae200db8`.
- Параметры: прямой Coulomb+LJ, вакуум, без периодики; выбран H40 как для всех трёх исходных поз.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-202/geometry-interaction.json`.
- Результат: OD2 — C4–O 3.375 Å, Hβ–O 2.573 Å, угол 129.22°; OD1 дальше. Прямой белок–лиганд Coulomb+LJ −140.64 кДж/моль.
- Отрицательный результат?: нет; прямой вакуумный вклад не учитывает растворитель и не является энергией связывания.
- Диагностика: модели содержат разное число воды; они не демонстрируют термодинамических популяций или сходимости.
- Влияние на следующие шаги: 202 имеет менее линейную геометрию, но численного критерия достаточности не установлено; нельзя считать позу исключённой только по углу.
- Статус: partial.

## 2026-11-01 13:32 | Минимизация комплекса позы 202
- Цель: исследовать вторую различимую позу в той же постановке, что кандидат 101.
- Скрипт/команда: `prepare_ksi_complex.py --protein enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb --ligand-sdf enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-202/ligand.sdf --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-202 --max-iterations 100 --padding-nm 0.8` в conda, CPU 2 потока.
- Входы: WT SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`, SDF `30b73607c4638da122c4a5aec0e024c62efda0655fe54132c70c6c6472da315a`.
- Параметры: Amber14/GAFF2.11/AM1-BCC, TIP3P, pH 7, ASH99 в обеих цепях, 100 итераций, 43 308 атомов.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-202/`, stdout `/tmp/opencode/ksi-complex-202.log`.
- Результат: общая потенциальная энергия −341435 → −619881 кДж/моль; размер системы отличается от кандидата 101 на 12 атомов.
- Отрицательный результат?: нет; общие энергии систем с разным количеством воды нельзя сравнивать для ранжирования связывания.
- Диагностика: вода сгенерирована заново, кристаллографические воды отсутствуют; сходимость и реальная структура не установлены.
- Влияние на следующие шаги: измерить одинаковые геометрические метрики; не делать вывод о сравнительной выгодности на основании общей энергии.
- Статус: partial.

## 2026-11-01 13:30 | Параметризация 5-AND для позы 202
- Цель: подготовить вторую различимую позу к физической релаксации.
- Скрипт/команда: `parameterize_ksi_ligand.py --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --xyz enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/candidate-202.xyz --output enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-202` в отдельном conda-окружении.
- Входы: 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`; XYZ `8c91372c9b301d1dc4b640553a00106f3a87632590f3a4d44c073bfa80dfb482`.
- Параметры: GAFF 2.11, AmberTools AM1-BCC, OpenMM CPU, 47 атомов, без растворителя.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-202/`.
- Результат: сумма зарядов ≈ 0; внутренняя энергия отдельного лиганда 187.66 → 162.65 кДж/моль.
- Отрицательный результат?: нет; это не энергия связывания и не полноценный комплекс.
- Диагностика: исходная поза 202 отличается от 101; зарядовый метод совпадает с предыдущей попыткой.
- Влияние на следующие шаги: можно построить отдельный комплекс той же физической модели.
- Статус: partial.

## 2026-11-01 13:29 | Причина пропущенного теста рендера
- Цель: подтвердить причину одного skipped из полного набора.
- Скрипт/команда: `python -m pytest -q -rs enzyme-kinetics/tests/test_render_ksi_preview.py`.
- Входы: `enzyme-kinetics/tests/test_render_ksi_preview.py`, проектное окружение Python 3.13.
- Параметры: однократный тест в core-окружении.
- Выходы: отчёт pytest.
- Результат: 1 skipped из-за отсутствия необязательного `imageio` в проектном окружении, pytest exit 5 поскольку не было выполненных тестов.
- Отрицательный результат?: да — тест рендера в core не выполнялся; реальный скрипт отдельно запущен в tools-окружении.
- Диагностика: отсутствие imageio, а не ошибка чтения файлов или графики.
- Влияние на следующие шаги: для теста установить pytest и запустить его в визуализационном окружении либо явно поддерживать отдельную группу зависимостей.
- Статус: partial.

## 2026-11-01 13:28 | Проверка нового диагностического кода
- Цель: выполнить полный набор тестов после добавления GPU, embedding, 3D-рендера и проверки атомов.
- Скрипт/команда: `python -m pytest -q` с проектным site-packages и OMP_NUM_THREADS=1.
- Входы: `enzyme-kinetics/tests/`, `enzyme-kinetics/scripts/`; индивидуальные хеши скриптов ещё не зафиксированы.
- Параметры: Python 3.13, CPU, ранее установленное проектное окружение.
- Выходы: вывод pytest.
- Результат: 64 passed, 1 skipped за 6.17 s.
- Отрицательный результат?: нет; пропущенный тест связан с отсутствующим imageio в проектном окружении; рендер отдельно проверен в tools-окружении.
- Диагностика: тесты проверяют программную логику, не доказывают качество GPU, белкового комплекса или механизма.
- Влияние на следующие шаги: регрессий не выявлено; научная валидация и GPU для полной системы остаются открытыми.
- Статус: done.

## 2026-11-01 13:27 | Проверка порядка атомов и стереохимии
- Цель: проверить перенос β-H через XYZ → SDF → минимизированный PDB.
- Скрипт/команда: `python enzyme-kinetics/scripts/verify_ksi_atom_mapping.py --sdf enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/ligand.sdf --xyz enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/candidate-101.xyz --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/minimized.pdb --assignment enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/atom-mapping-check.json`.
- Входы: SDF SHA-256 `752c8bcbf25d6fd60bbde6696d0183fc9de647e89f2970e4304fd97844725b19`; XYZ `d454619a37b907ac61150f95a38f753b82c9b494b8cd9f2e9d67e35fc01cd862`; PDB `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`; assignment `794d5149e71874cee9093a1af456559eb7c8c33abdea4a7b6f150a28ae200db8`.
- Параметры: RDKit 2026.03.1, OpenMM PDB parser; пересчёт стереохимии из 3D и проверка связи по SDF.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/atom-mapping-check.json`.
- Результат: максимальная разница XYZ–SDF 0.0000704 Å; атом 40 H связан в SDF с C4 индекс 13; C4–H 1.0963 Å до/после; изомерный SMILES из 3D сохраняется после минимизации.
- Отрицательный результат?: нет; тождество одинаковых атомов в PDB всё ещё опирается на сохранённый порядок writer, изотопно/экспериментально перенос не измерен.
- Диагностика: проверка связи и 3D-стереохимии независима от предыдущего сравнения только символов XYZ; связь реальной позы с экспериментом не доказывается.
- Влияние на следующие шаги: повышена уверенность в модельном индексе β-H; дальнейшие QM/MM требуют дополнительной проверки структуры и ансамбля.
- Статус: partial для реального механизма, done для проверенных файлов.

## 2026-11-01 13:25 | 3D-визуализация модельного комплекса
- Цель: проверить молекулярный рендер и экспорт MP4 на реальных координатах.
- Скрипт/команда: `python enzyme-kinetics/scripts/render_ksi_preview.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/minimized.pdb --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/visualization-preview --frames 24`.
- Входы: PDB SHA-256 `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`.
- Параметры: matplotlib Agg, Bio.PDB, imageio/FFmpeg, 24 кадра 8 fps; вода отфильтрована только для рендера, координаты не менялись.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/visualization-preview/preview.png`, `camera-orbit.mp4`, `summary.json`.
- Результат: визуализированы 125 CA каждой цепи и 21 тяжёлый атом стероида, MP4 создан; FFmpeg масштабировал кадры 600×600 до 608×608.
- Отрицательный результат?: нет; это вращение камеры вокруг одной статичной модельной структуры, не MD и не видео химической реакции.
- Диагностика: решена проблема hybrid-36 путём исключения воды и ионов перед Bio.PDB; исходный файл не изменён.
- Влияние на следующие шаги: функциональный 3D-рендер доступен; реалистичная анимация возможна лишь из валидированной траектории.
- Статус: done для визуализации, не для динамики реакции.

## 2026-11-01 13:23 | Неудачная попытка 3D-рендера
- Цель: визуализировать полученный молекулярный комплекс.
- Скрипт/команда: `python enzyme-kinetics/scripts/render_ksi_preview.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/minimized.pdb --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/visualization-preview --frames 24`.
- Входы: PDB SHA-256 `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`.
- Параметры: Bio.PDB, matplotlib Agg, imageio/FFmpeg; 24 кадра статичной модели.
- Выходы: traceback, папка не создана.
- Результат: Bio.PDB не разобрал hybrid-36 residue ID `A000` у большого количества растворителя.
- Отрицательный результат?: да — MP4 и PNG не созданы.
- Диагностика: рендер не требует воды; безопасно отфильтровать только записи белка и лиганда до разбора, не меняя оригинальный PDB.
- Влияние на следующие шаги: исправить чтение и повторить в новом запуске.
- Статус: blocked до исправления.

## 2026-11-01 13:21 | PySCF: проверка электростатического embedding
- Цель: проверить квантовый расчёт с фиксированным внешним зарядом.
- Скрипт/команда: `python enzyme-kinetics/scripts/check_qmmm_embedding.py --output enzyme-kinetics/outputs/active-site-v3/qmmm-electrostatic-smoke.json`.
- Входы: синтетическая молекула H2O и точечный заряд +0.5e, координаты заданы в скрипте; внешние файлы отсутствуют.
- Параметры: PySCF 2.14.0, RHF/STO-3G, CPU 1 поток, заряд в (3,0,0) Å.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/qmmm-electrostatic-smoke.json`.
- Результат: оба SCF сошлись; разница электронных энергий −0.00273552 Hartree.
- Отрицательный результат?: нет для теста; это не связанный PySCF/OpenMM QM/MM драйвер и не ферментативная реакция.
- Диагностика: проверяется только электростатическое влияние одного внешнего заряда, без QM/MM-границы, сил MM, поляризации и динамики.
- Влияние на следующие шаги: базовый компонент доступен, но полноценное QM/MM и путь реакции не подтверждены.
- Статус: partial.

## 2026-11-01 13:19 | Полный комплекс: CPU PME
- Цель: проверить, устраняет ли перенос PME на CPU ошибку GPU для полного комплекса.
- Скрипт/команда: `run-rusticl-openmm.sh validate_opencl_gpu.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/system.xml --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/opencl-validation-cpupme.json --threads 2 --steps 2 --precision single --cpu-pme`.
- Входы: PDB SHA-256 `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`, XML `7853e0ae0a3a43067a59e2fc3b9f6dd9e94121f3d3d6411fb49eb29a005517a4`.
- Параметры: OpenMM 8.6.1, OpenCL single, `UseCpuPme=true`, два потока CPU.
- Выходы: JSON не создан.
- Результат: предупреждение неподдерживаемой OpenCL-реализации и аварийный exit 139.
- Отрицательный результат?: да — отключение GPU PME не устранило аварию; истинная причина пока не установлена.
- Диагностика: тот же лигандный тест без периодической среды работал; это не доказывает, что сам PME является причиной.
- Влияние на следующие шаги: полный комплекс считать несовместимым с текущей комбинацией Rusticl и OpenMM до дальнейшей диагностики.
- Статус: blocked для расчёта полного комплекса на GPU.

## 2026-11-01 13:17 | CPU/GPU: изолированный 5-AND
- Цель: отделить проблему полного комплекса от работоспособности малого реального расчёта.
- Скрипт/команда: `bash enzyme-kinetics/scripts/run-rusticl-openmm.sh enzyme-kinetics/scripts/validate_opencl_gpu.py --pdb enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/ligand-minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/ligand-system.xml --output enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/opencl-single.json --threads 2 --steps 2 --precision single`.
- Входы: PDB SHA-256 `000782c1fcbf7601c2a8449a3e29a6f9c0019ece5ce220fd2d0d5a8bf5d0063d`, XML `8e1b4742b79635a263495aba075a7ac880c1fc04bfb0d3bc4decc15fd576929e`.
- Параметры: OpenMM 8.6.1, CPU 2 потока, Rusticl single, 47 атомов, два шага.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/opencl-single.json`.
- Результат: CPU/GPU ΔE −0.0000263 кДж/моль, RMS-компонент силы 0.04255 и максимум 0.1973 кДж/(моль·нм); GPU выполнил два шага.
- Отрицательный результат?: нет для малого теста; OpenMM предупреждает об unsupported реализации, полный комплекс ранее аварийно завершился.
- Диагностика: сравнение сделано по одной конфигурации без растворителя и PME; результат не подтверждает длительную устойчивость или точность всех силовых компонентов.
- Влияние на следующие шаги: GPU работает для малого случая, но к полной задаче применять только после локализации аварии.
- Статус: partial.

## 2026-11-01 13:15 | Single и mixed на Radeon
- Цель: проверить режимы точности OpenCL на устройстве и полном комплексе.
- Скрипт/команда: минимальный OpenMM Context с `Precision=mixed`; затем `run-rusticl-openmm.sh validate_opencl_gpu.py` для 43 320 атомов с `--precision single --threads 2 --steps 2`.
- Входы: PDB SHA-256 `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`, XML SHA-256 `7853e0ae0a3a43067a59e2fc3b9f6dd9e94121f3d3d6411fb49eb29a005517a4`; исправленные локальные headers libclang 20.
- Параметры: OpenMM 8.6.1, Mesa Rusticl 25.2.8, Radeon RX 550, isolated bwrap.
- Выходы: JSON не создан; вывод и код возврата.
- Результат: простой mixed контекст не создался (`No compatible OpenCL platform`); единичное простое ядро single ранее вернуло 1.0 кДж/моль; полный комплекс в single аварийно завершился с кодом 139.
- Отрицательный результат?: да — реальная система на GPU сейчас нестабильна; OpenMM предупреждает, что реализация не поддерживается.
- Диагностика: смешанная точность недоступна; причина segfault ещё не локализована. Изолированная пользовательская установка не изменила системные драйверы.
- Влияние на следующие шаги: не запускать production и QM/MM на GPU до выявления причины и сравнения с CPU.
- Статус: blocked для полного комплекса.

## 2026-11-01 13:12 | Пробная проверка OpenCL на полном комплексе
- Цель: сравнить силы и энергию OpenMM CPU и RX 550 для 43 320-атомного комплекса.
- Скрипт/команда: `bash enzyme-kinetics/scripts/run-rusticl-openmm.sh enzyme-kinetics/scripts/validate_opencl_gpu.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/system.xml --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/opencl-validation.json --threads 2 --steps 2`.
- Входы: PDB SHA-256 `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`; system XML SHA-256 `7853e0ae0a3a43067a59e2fc3b9f6dd9e94121f3d3d6411fb49eb29a005517a4`.
- Параметры: OpenMM 8.6.1, CPU 2 потока, OpenCL mixed, два шага; изолированный Rusticl.
- Выходы: traceback; JSON не создан.
- Результат: CPU этап выполнен, создание второго контекста OpenCL завершилось `No compatible OpenCL platform is available`, хотя в отдельном процессе простое ядро запускалось.
- Отрицательный результат?: да — CPU/GPU сравнение пока не получено.
- Диагностика: проверить отличия wrapper и запуск в отдельных процессах, ресурсные инициализации платформ; не считать GPU корректным на полном комплексе.
- Влияние на следующие шаги: GPU остаётся экспериментальным и не используется для научных выводов.
- Статус: partial.

## 2026-11-01 13:08 | Rusticl обнаружен, ядро не скомпилировалось
- Цель: запустить OpenMM на Radeon RX 550 без изменения системных файлов.
- Скрипт/команда: локальная распаковка Ubuntu-пакетов mesa-opencl-icd/libclc-20/libclang-cpp20/libllvmspirvlib20.1/clinfo; `bwrap` с пользовательским ICD, `RUSTICL_ENABLE=radeonsi` и монтированием libclc; затем создание простого OpenCL Context.
- Входы: Mesa пакет SHA-256 `07d3d2189823918c005a24b9d2c73e0bed0ebee70cd08dd7584a4f649247067c`; libclc-20 SHA-256 `588595000d6b335120ca71c22a9d366b302c17fde2573f525b844846b374924b`; остальные пакеты Ubuntu репозитория.
- Параметры: Mesa Rusticl 25.2.8, amdgpu, RX 550 polaris12, OpenMM 8.6.1; bwrap без sudo.
- Выходы: пользовательский каталог `$HOME/.local/share/ksi-mesa-opencl/`; вывод clinfo и traceback.
- Результат: clinfo обнаружил GPU; OpenMM видит OpenCL, но создание Context закончилось `opencl-c-base.h file not found`; предупреждение OpenMM: unsupported OpenCL implementation.
- Отрицательный результат?: да — реальный расчёт пока не выполнен. Первая попытка bwrap с read-only `/usr/lib` также не смогла создать точку монтирования; повтор с отдельным tmpfs показал устройство.
- Диагностика: libclc виден после bind-mount, но clang resource headers отсутствуют; необходимо установить их локально и повторить.
- Влияние на следующие шаги: определить совместимость OpenMM после устранения отсутствующего заголовка; не объявлять обнаружение GPU успешным вычислением.
- Статус: partial.

## 2026-11-01 11:52 | GPU и OpenCL: первичная диагностика
- Цель: установить причину отсутствия GPU в OpenMM без изменений системы.
- Скрипт/команда: `lspci -nnk -s 03:00.0`, `id`, `getfacl /dev/dri/renderD128 /dev/kfd`, `apt-cache policy mesa-opencl-icd`, `openmm.Platform.getPluginLoadFailures()`, проверка `/etc/OpenCL/vendors` и `os.access`.
- Входы: PCI и системные устройства, OpenMM из `$HOME/.local/share/enz-kinetics-mm/`; хеши изменяемого состояния системы не применимы.
- Параметры: Linux 6.17, OpenMM 8.6.1, пользовательский процесс без sudo.
- Выходы: журнал; системные настройки не менялись.
- Результат: AMD Lexa PRO `[1002:699f]`, драйвер amdgpu; renderD128 доступен по ACL, /dev/kfd недоступен; OpenCL vendor ICD отсутствует, Mesa OpenCL пакет не установлен; OpenMM OpenCL plugin присутствует, но платформы только CPU/Reference.
- Отрицательный результат?: да — GPU расчёт пока невозможен; NVIDIA CUDA отсутствует ожидаемо, HIP plugin не загружается из-за отсутствия libhiprtc.
- Диагностика: причина отсутствия OpenCL конкретизирована до отсутствующего ICD, доступности /dev/kfd и неизвестной совместимости GPU/runtime; системные исправления пока не выполнялись.
- Влияние на следующие шаги: безопасно испытать локальную Mesa OpenCL, не менять группы/драйверы без согласования.
- Статус: partial. Время по системным часам; ранее некоторые записи имеют более позднее время.

## 2026-11-01 11:56 | Проверка кодирования анимации
- Цель: проверить локальное создание MP4 без GPU.
- Скрипт/команда: Python imageio `get_writer(..., fps=5, codec="libx264")`, 8 синтетических RGB-кадров 64×64.
- Входы: синтетические кадры numpy, внешние файлы отсутствуют.
- Параметры: imageio 2.38.1, imageio-ffmpeg 0.7.0, FFmpeg 7.0.2-static, libx264.
- Выходы: `/tmp/opencode/ksi-animation-check.mp4`.
- Результат: writer завершился без ошибки.
- Отрицательный результат?: нет; молекулярный 3D-рендер, видео из траектории и воспроизведение MP4 отдельно не проверены.
- Диагностика: тестирует только программный кодировщик и простые кадры.
- Влияние на следующие шаги: для научной анимации ещё нужны корректные траектории и проверка 3D-рендерера.
- Статус: done для кодирования, partial для полноценной визуализации.

## 2026-11-01 11:55 | Экспорт версий инструментов
- Цель: сохранить состав двух изолированных окружений.
- Скрипт/команда: попытка `python -m pip freeze`, затем `uv pip freeze --python <tools-python>`; `micromamba list -p <mm-prefix> --explicit`.
- Входы: окружения `$HOME/.local/share/enzyme-kinetics-tools/` и `$HOME/.local/share/enz-kinetics-mm/`.
- Параметры: pip-окружение Python 3.12 и conda-forge Python 3.11.
- Выходы: `docs/experiments/active-site-v3/tools-pip-freeze.txt`, `docs/experiments/active-site-v3/tools-conda-explicit.txt`.
- Результат: сохранены версии Python-пакетов и явные conda URLs.
- Отрицательный результат?: да — первый вызов `python -m pip` не сработал, поскольку uv venv не содержит pip; `uv pip freeze` успешен.
- Диагностика: замена команды не меняла установленные пакеты; conda explicit list не гарантирует вечную доступность URL.
- Влияние на следующие шаги: зафиксированы составы окружений; проектный uv.lock не менялся.
- Статус: done с зафиксированной неудачной попыткой.

## 2026-11-01 11:54 | Тесты β-H и параметризации
- Цель: проверить новый код и регрессии проекта.
- Скрипт/команда: `python -m pytest -q`.
- Входы: скрипты `assign_ksi_beta_h.py` SHA-256 `389d90ac7a598154622fd56c3f78b36721aaa60172a447a071d9ffab9bcf87d9`, `parameterize_ksi_ligand.py` `c68034f69c624fe0366242626fc0014289b7bbf02669730b350d1d6f29fba755`, `prepare_ksi_complex.py` `6bb3f599faecf87873513aadbcd3432e4aa1cdaf71f03db76dcfc2880d289351`, `analyze_ksi_complex.py` `5a0272e87c4db5ae66bd872e134a8d02435897dbc0667e9af0f33e61d24ab69e` в `enzyme-kinetics/scripts/`.
- Параметры: Python 3.13 с ранее установленными site-packages; функциональный расчёт OpenFF выполнялся отдельно в conda.
- Выходы: вывод pytest.
- Результат: 61 passed за 5.83 s.
- Отрицательный результат?: нет; unit-тесты не валидируют научную модель комплекса или точность энергии связывания.
- Диагностика: проверены ориентация β-H при отражении, вырожденная геометрия, индексирование XYZ, варианты ASH и формула прямого взаимодействия.
- Влияние на следующие шаги: конвейер проверен программно; необходима физическая и экспериментальная валидация.
- Статус: done.

## 2026-11-01 11:52 | Геометрия и прямые взаимодействия комплекса
- Цель: проверить модельный β-H после минимизации и приблизительную прямую энергию белок–лиганд.
- Скрипт/команда: `python enzyme-kinetics/scripts/analyze_ksi_complex.py --pdb enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/minimized.pdb --system enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/system.xml --beta-assignment enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/geometry-interaction.json`.
- Входы: PDB SHA-256 `fd59bb580ab9ae0a689daa8f776224c8e23acc708a3c35ee17676bd01a0befd8`; XML `7853e0ae0a3a43067a59e2fc3b9f6dd9e94121f3d3d6411fb49eb29a005517a4`; β-H JSON `794d5149e71874cee9093a1af456559eb7c8c33abdea4a7b6f150a28ae200db8`.
- Параметры: парные Coulomb + Lennard-Jones, без периодических образов, без отсечения и без экранирования растворителем; силы и заряды из системы Amber14/GAFF.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/geometry-interaction.json`.
- Результат: для OD1 C4–O 3.530 Å, Hβ–O 2.470 Å, C4–Hβ–O 162.11°; прямой парный Coulomb+LJ белок–лиганд −141.36 кДж/моль.
- Отрицательный результат?: нет, но указанная величина не PME-энергия, не учитывает экранирование, десольватацию и энтропию и **не является энергией связывания**; одной позы недостаточно для вывода.
- Диагностика: 47 атомов лиганда, 3760 белка; полная система содержит воду, но анализ прямой пары её исключает; β-H присваивается по сохранённому порядку атомов, требуется независимая графовая проверка.
- Влияние на следующие шаги: геометрия после короткой минимизации совместима с предварительным фильтром; нужны альтернативы и корректная оценка неопределённости.
- Статус: partial — исследовательская оценка, не валидация комплекса.

## 2026-11-01 11:49 | Минимизация WT–5-AND в TIP3P
- Цель: пробная физическая параметризация комплекса и короткая релаксация кандидата 101.
- Скрипт/команда: `python enzyme-kinetics/scripts/prepare_ksi_complex.py --protein enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb --ligand-sdf enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/ligand.sdf --output enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101 --max-iterations 100 --padding-nm 0.8`, PATH conda, `OPENMM_CPU_THREADS=2`.
- Входы: WT SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`; ligand.sdf SHA-256 `752c8bcbf25d6fd60bbde6696d0183fc9de647e89f2970e4304fd97844725b19`.
- Параметры: Amber14 белок, GAFF 2.11/AM1-BCC лиганд, TIP3P, PME 0.8 nm, принудительный ASH99 в обеих цепях как гипотеза, остальные варианты OpenMM pH 7, нейтрализация, 43 320 атомов, 100 итераций; кристаллографические воды не сохранены.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/complex-gaff2-101/` — `minimized.pdb`, `system.xml`, `summary.json`.
- Результат: потенциальная энергия системы −343942 → −617061 кДж/моль; неприсутствующих остатков и тяжёлых атомов не достраивали.
- Отрицательный результат?: нет, но общая потенциальная энергия растворённой системы не является энергией взаимодействия или связывания; сходимость и валидность позы не установлены.
- Диагностика: SDF и WT входы захешированы; ограниченная минимизация без MD, только одна модель протонирования и один кандидат. Равновесность и сохранение β-H геометрии требуют проверки.
- Влияние на следующие шаги: сравнить координаты до/после, идентичность лиганда и β-H; рассчитать корректно определённую взаимодействующую энергию только с явными оговорками и согласованными условиями.
- Статус: partial — пробная физическая модель, не экспериментально подтверждённый WT–ES и не энергия связывания.

## 2026-11-01 11:46 | GAFF2/AM1-BCC: изолированный 5-AND
- Цель: проверить физическую параметризацию и минимизацию отдельного лиганда seed 101.
- Скрипт/команда: `python enzyme-kinetics/scripts/parameterize_ksi_ligand.py --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --xyz enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/candidate-101.xyz --output enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2` с PATH conda и `OPENMM_CPU_THREADS=2`.
- Входы: 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`; XYZ SHA-256 `d454619a37b907ac61150f95a38f753b82c9b494b8cd9f2e9d67e35fc01cd862`.
- Параметры: AmberTools AM1-BCC явно выбран, GAFF 2.11, OpenMM 8.6.1 CPU, 200 итераций, NoCutoff, без белка и растворителя.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101-v2/` — SDF, минимизированный PDB, OpenMM system XML, summary.
- Результат: 47 атомов, сумма частичных зарядов ≈ 1.3e-16 e; внутренняя потенциальная энергия 187.66 → 162.66 кДж/моль.
- Отрицательный результат?: нет; энергия не является энергией взаимодействия и не относится к комплексу.
- Диагностика: исправлены PDB residue IDs и выбор toolkit после неуспешной попытки 11:44; геометрия ES и порядок атомов требуют отдельной проверки.
- Влияние на следующие шаги: backend лигандной параметризации функционирует, но сборка и краткая минимизация белок–лиганд–вода ещё необходимы.
- Статус: partial — только изолированный лиганд, а не WT–ES.

## 2026-11-01 11:44 | Первая попытка GAFF2/AM1-BCC 5-AND
- Цель: получить параметры и минимизировать отдельный лиганд seed 101.
- Скрипт/команда: `python enzyme-kinetics/scripts/parameterize_ksi_ligand.py --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --xyz enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/candidate-101.xyz --output enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101` с PATH conda и `OPENMM_CPU_THREADS=2`.
- Входы: 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`; XYZ SHA-256 `d454619a37b907ac61150f95a38f753b82c9b494b8cd9f2e9d67e35fc01cd862`.
- Параметры: OpenFF 0.18.0, AmberTools, GAFF 2.11, OpenMM CPU; ligand-only.
- Выходы: неполный каталог `enzyme-kinetics/outputs/active-site-v3/ligand-gaff2-am1bcc-101/`.
- Результат: создание и минимизация системы прошли до записи, однако PDB writer вызвал TypeError из-за числового residue id; JSON и XML отсутствуют.
- Отрицательный результат?: да — попытка завершилась исключением; предупреждение NAGL также показывает, что registry пытался другие backend, так что выбор AmberTools следует закрепить явно.
- Диагностика: ошибка в `PDBFile.writeFile(..., keepIds=True)`; использовать автогенерируемые ID и явный toolkit wrapper.
- Влияние на следующие шаги: сохранить неполные файлы и повторить в новом каталоге с исправлениями.
- Статус: blocked до повтора.

## 2026-11-01 11:41 | Проверка backend AmberTools
- Цель: проверить доступность инструмента AM1-BCC, а не только факт установки.
- Скрипт/команда: импорт OpenFF, OpenMM, openmmforcefields, pdbfixer; `AmberToolsToolkitWrapper.is_available()` сначала без, затем с добавлением `enz-kinetics-mm/bin` в PATH процесса.
- Входы: окружение `$HOME/.local/share/enz-kinetics-mm/`; версии каталогов через conda-forge.
- Параметры: OpenFF 0.18.0, OpenMM 8.6.1, openmmforcefields 0.16.0, Python 3.11.
- Выходы: вывод команд; настройки оболочки не менялись.
- Результат: backend без PATH не найден; с PATH доступен (`True`), модули импортируются.
- Отрицательный результат?: да — автоматическая доступность без PATH отсутствует; до успешного расчёта заряды и параметры не подтверждены.
- Диагностика: antechamber установлен внутри conda-окружения и требует добавления его bin к PATH.
- Влияние на следующие шаги: запускать расчёты с явным PATH и проверить AM1-BCC на конкретном 5-AND.
- Статус: partial до теста реального лиганда.

## 2026-11-01 11:40 | Повтор установки conda-forge
- Цель: завершить прерванную установку физической параметризации.
- Скрипт/команда: `micromamba install -y -p "$HOME/.local/share/enz-kinetics-mm" -c conda-forge --override-channels python=3.11 openff-toolkit openmmforcefields ambertools openmm rdkit pdbfixer`.
- Входы: conda-forge repodata и ранее загруженный cache; micromamba SHA-256 `701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd`.
- Параметры: отдельный MAMBA_ROOT_PREFIX, CPU, таймаут 600 s, без изменения проектного окружения.
- Выходы: `$HOME/.local/share/enz-kinetics-mm/`, текстовый лог `/tmp/opencode/ksi-mamba-install.log`.
- Результат: команда вернула код 0; окончательные версии и работоспособность проверить отдельно.
- Отрицательный результат?: нет; установка не доказывает параметризацию или физическую корректность комплекса.
- Диагностика: повтор транзакции вместо удаления частично загруженных файлов.
- Влияние на следующие шаги: проверка backend и параметризация отдельного лиганда.
- Статус: partial до функционального теста.

## 2026-11-01 11:34 | Прерванная установка conda-forge
- Цель: установить проверенный стек OpenFF/AmberTools для параметризации лигандов.
- Скрипт/команда: `micromamba create -y -p "$HOME/.local/share/enz-kinetics-mm" -c conda-forge --override-channels python=3.11 openff-toolkit openmmforcefields ambertools openmm rdkit pdbfixer` с отдельным MAMBA_ROOT_PREFIX.
- Входы: официальные conda-forge repodata; micromamba SHA-256 `701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd`.
- Параметры: CPU, Python 3.11, план 248 пакетов и 936 MB загрузки; проектная .venv не менялась.
- Выходы: возможные неполные каталоги `$HOME/.local/share/enz-kinetics-mm/` и `$HOME/.local/share/enz-kinetics-mamba/`.
- Результат: зависимости разрешены; установка прервана общим таймаутом 120 s после начала транзакции.
- Отрицательный результат?: да — завершённость и работоспособность окружения не подтверждены.
- Диагностика: необходимо проверить состояние, затем безопасно повторить транзакцию без удаления исходных файлов.
- Влияние на следующие шаги: до успешной проверки не считать параметризацию доступной.
- Статус: partial, причина — таймаут.

## 2026-11-01 11:31 | Загрузка micromamba
- Цель: получить изолированный менеджер для официального conda-forge backend OpenFF/AmberTools.
- Скрипт/команда: `curl -L --fail https://micro.mamba.pm/api/micromamba/linux-64/latest -o /tmp/opencode/micromamba-package.tar.bz2`; `tar -xjf ... -C /tmp/opencode bin/micromamba`.
- Входы: официальный endpoint micromamba; архив SHA-256 `9a97976244929bbae28fb1fae504e74a3e7e6d1932870e3a03321f01b7fbd07b`.
- Параметры: Linux x86-64; без изменения настроек оболочки и проектной `.venv`.
- Выходы: `/tmp/opencode/bin/micromamba`, SHA-256 `701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd`.
- Результат: бинарник извлечён; официальная независимая сверка checksum UNKNOWN.
- Отрицательный результат?: нет; до создания окружения OpenFF не установлен.
- Диагностика: использован официальный документированный URL, содержимое только временного каталога.
- Влияние на следующие шаги: создание отдельного conda-forge окружения, если разрешение зависимостей и ресурсный бюджет позволят.
- Статус: partial.

## 2026-11-01 11:29 | Геометрическое назначение β-H
- Цель: различить два H при C4 в трёх кандидатах.
- Скрипт/команда: `python enzyme-kinetics/scripts/assign_ksi_beta_h.py --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --pose-dir enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water --output enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json`.
- Входы: 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`; три XYZ и их SHA-256 в выходном JSON и записи 11:14.
- Параметры: RDKit 2026.09.1; сторона метильной группы C19 определена как β, локальная плоскость C3/C5/C10, использованы координаты кандидатов.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/ksi-beta-h-assignment.json`.
- Результат: для всех трёх кандидатов H с индексом 40 на стороне C19 (геометрически β), H41 — на противоположной; индексы только для этих XYZ.
- Отрицательный результат?: нет; это модельное назначение не доказывает экспериментальную геометрию комплекса, XYZ не содержит связей и не гарантирует корректный порядок одних лишь одинаковых элементов.
- Диагностика: знаки проекций C10→C19 и C4→H40 согласованы, проекция H41 противоположна; отдельное графовое/визуальное подтверждение требуется до QM/MM.
- Влияние на следующие шаги: сопоставить β-H с геометрией O(Asp38), но не использовать назначение как доказательство механизма или комплекса.
- Статус: partial до независимой проверки атомного порядка и стереохимии.

## 2026-11-01 11:27 | Установка openmmforcefields отдельно
- Цель: проверить наличие отдельного пакета генераторов силовых полей.
- Скрипт/команда: `uv pip install --python "$HOME/.local/share/enzyme-kinetics-tools/bin/python" openmmforcefields`.
- Входы: PyPI; хеш дистрибутива UNKNOWN.
- Параметры: Python 3.12, установка без openff-toolkit.
- Выходы: пользовательское окружение расчётных инструментов.
- Результат: установлен openmmforcefields 0.15.1.
- Отрицательный результат?: нет; установка не включает автоматически OpenFF Toolkit или AmberTools и не означает доступности полного backend.
- Диагностика: предыдущий общий запрос был заблокирован openff-toolkit; отдельная установка пакета прошла.
- Влияние на следующие шаги: необходимо проверять реальную параметризацию стероида, а не только импорт.
- Статус: partial.

## 2026-11-01 11:25 | Установка OpenFF через PyPI не удалась
- Цель: установить параметризацию лиганда для физической минимизации комплекса.
- Скрипт/команда: `uv pip install --python "$HOME/.local/share/enzyme-kinetics-tools/bin/python" openmmforcefields openff-toolkit`.
- Входы: PyPI package index; локальный пользовательский Python 3.12.
- Параметры: общий запрос двух пакетов, CPU.
- Выходы: traceback резольвера; установка не произведена.
- Результат: резольвер сообщил, что единственная найденная версия openff-toolkit на PyPI отозвана; набор зависимостей неразрешим.
- Отрицательный результат?: да — OpenFF/GAFF параметризация в этом окружении не установлена.
- Диагностика: в сообщении резольвера причина yanked — некоординированный релиз на PyPI; не использовать отозванную сборку.
- Влияние на следующие шаги: требуется официальный conda-forge или альтернативный проверенный путь; произвольно назначать параметры нельзя.
- Статус: blocked до доступного силового поля и charge backend.

## 2026-11-01 11:24 | Исправление зависимости Meeko
- Цель: устранить ошибку импорта и повторить smoke-тест.
- Скрипт/команда: `uv pip install --python "$HOME/.local/share/enzyme-kinetics-tools/bin/python" gemmi`; затем Python-импорты Vina, Meeko, RDKit, OpenMM, ASE, PySCF, tblite, py3Dmol, imageio и imageio_ffmpeg.
- Входы: PyPI gemmi и установленное пользовательское окружение; хеши дистрибутивов UNKNOWN.
- Параметры: gemmi 0.7.5, Python 3.12.
- Выходы: обновлённое окружение.
- Результат: импорты успешны, FFmpeg 7.0.2-static, OpenMM platforms Reference и CPU.
- Отрицательный результат?: нет; GPU-платформа отсутствует, полноценный расчёт и рендер не проверены.
- Диагностика: повторён тот же тест, который завершился на отсутствии gemmi.
- Влияние на следующие шаги: подготовка и минимальные функциональные проверки.
- Статус: done для импорта, partial для полной функциональности.

## 2026-11-01 11:23 | Неудачный smoke-тест Meeko
- Цель: проверить импорты установленных инструментов и наличие FFmpeg.
- Скрипт/команда: Python-импорт `vina,meeko,rdkit,openmm,ase,pyscf,tblite,py3Dmol,imageio,imageio_ffmpeg`.
- Входы: установленное окружение `$HOME/.local/share/enzyme-kinetics-tools/`; checksum дистрибутивов UNKNOWN.
- Параметры: Python 3.12, локальный запуск.
- Выходы: traceback в журнале текущей сессии.
- Результат: импорт остановился на `ModuleNotFoundError: No module named gemmi` внутри Meeko; дальнейшие импорты этим запуском не проверены.
- Отрицательный результат?: да — одного пакета Meeko недостаточно для импорта.
- Диагностика: трассировка указывает на `meeko/chemtempgen.py` и зависимость gemmi.
- Влияние на следующие шаги: установить недостающий gemmi и повторить проверку.
- Статус: blocked до исправления зависимости.

## 2026-11-01 11:22 | Установка вычислительных инструментов
- Цель: локальные инструменты докинга, QM и визуализации без вмешательства в проектный lock.
- Скрипт/команда: `uv pip install --python "$HOME/.local/share/enzyme-kinetics-tools/bin/python" vina meeko rdkit openmm ase pyscf tblite py3Dmol imageio imageio-ffmpeg biopython scipy matplotlib`.
- Входы: PyPI-пакеты; lock до установки отсутствовал, исходные checksum пакетов UNKNOWN.
- Параметры: Python 3.12 в отдельном пользовательском окружении, CPU.
- Выходы: `$HOME/.local/share/enzyme-kinetics-tools/`.
- Результат: установлены Vina 1.2.7, Meeko 0.8.0, RDKit 2026.9.1, OpenMM 8.6.1, PySCF 2.14.0, tblite 0.8.0, ASE 3.29.0, py3Dmol 2.5.5, imageio 2.38.1, imageio-ffmpeg 0.7.0 и зависимости.
- Отрицательный результат?: нет; функциональность и комплекты силовых полей ещё требуют проверки; готовность к QM/MM не следует из установки пакетов.
- Диагностика: пакеты установлены вне `.venv`; сетевой источник PyPI, проверка воспроизводимости будет зафиксирована отдельно.
- Влияние на следующие шаги: smoke-тест импорта, визуализации и подготовки лигандов.
- Статус: partial до smoke-теста.

## 2026-11-01 11:21 | Изолированное окружение расчётных инструментов
- Цель: установить инструменты, не изменяя окружение и lock проекта.
- Скрипт/команда: `uv venv --python /usr/bin/python3.12 "$HOME/.local/share/enzyme-kinetics-tools"`.
- Входы: системный Python 3.12.3, pyproject.toml; бинарный SHA-256 UNKNOWN.
- Параметры: отдельное пользовательское виртуальное окружение, вне `.venv` проекта.
- Выходы: `$HOME/.local/share/enzyme-kinetics-tools/`.
- Результат: окружение создано; uv предупреждает о несовместимости Python 3.12 с требованием 3.13 проекта, поэтому оно предназначено для отдельных инструментов, не для запуска самого проекта.
- Отрицательный результат?: нет; воспроизводимость зависит от зафиксированных версий установленных пакетов.
- Диагностика: действующее окружение `.venv` не перезаписывалось.
- Влияние на следующие шаги: отдельная установка докинга и визуализации.
- Статус: done.

## 2026-11-01 11:15 | Тесты геометрического предскрининга
- Цель: проверить проект и вычисление геометрии после добавления скрипта.
- Скрипт/команда: `python -m pytest -q`.
- Входы: `enzyme-kinetics/scripts/screen_ksi_poses.py` SHA-256 `ef399f2d599508269a4ca364ba21b139fa8d2d72a8217f13ec36c9d8566eb0f9`; `enzyme-kinetics/tests/test_screen_ksi_poses.py` SHA-256 `513b72731defb1ac1481600a3a101cabb933819f42103e395d30ee47f3a6e4c8`.
- Параметры: Python 3.13 и существующие site-packages.
- Выходы: вывод pytest.
- Результат: 56 passed за 5.57 s.
- Отрицательный результат?: нет.
- Диагностика: тест геометрии проверяет угол и расстояния, но не проверяет полноту поиска или физическую достоверность комплекса.
- Влияние на следующие шаги: программный контроль прошёл, научная валидация остаётся открытой.
- Статус: done.

## 2026-11-01 11:14 | Дополнение к записи от 11:13: SHA-256
- Цель: зафиксировать хеши входов сравнения поз.
- Скрипт/команда: `sha256sum enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/candidate-{101,202,303}.xyz`.
- Входы: `candidate-101.xyz` `d454619a37b907ac61150f95a38f753b82c9b494b8cd9f2e9d67e35fc01cd862`; `candidate-202.xyz` `8c91372c9b301d1dc4b640553a00106f3a87632590f3a4d44c073bfa80dfb482`; `candidate-303.xyz` `23f82d34086c3b246aec49fbfb53c11c4eb0b4e3bd258cc1ec234d17bba1995c` — все пути относительно `enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/`.
- Параметры: SHA-256 файлов XYZ.
- Выходы: эта запись журнала.
- Результат: хеши зафиксированы; отдельного файла хешей в каталоге не создавалось.
- Отрицательный результат?: нет.
- Диагностика: исправление указания предыдущей записи на якобы имеющийся локальный каталог хешей.
- Влияние на следующие шаги: обеспечено воспроизведение сравнения RMSD.
- Статус: done.

## 2026-11-01 11:13 | Различие поз с водой
- Цель: выяснить, получили ли три старта различимые ориентации.
- Скрипт/команда: Python stdlib, попарный RMSD 21 тяжёлого атома из `candidate-{101,202,303}.xyz` в неизменённой системе координат WT.
- Входы: `enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/candidate-{101,202,303}.xyz`, хеши приведены в локальном каталоге с результатом; случайные seed 101, 202, 303.
- Параметры: сравнение без отдельного выравнивания лиганда, что сохраняет различия ориентации в активном центре.
- Выходы: вывод команды.
- Результат: 101–303: 0.000106 Å; 101–202 и 303–202: 2.432 Å.
- Отрицательный результат?: да — три старта дали только два различных геометрических минимума, не три независимых состояния комплекса.
- Диагностика: совпадение 101/303 превышает точность записи XYZ; разнообразие двух оставшихся не указывает на термодинамическую популяцию.
- Влияние на следующие шаги: не объявлять набор подтверждённым ансамблем; требуются новые начальные гипотезы и физическая релаксация.
- Статус: done.

## 2026-11-01 11:12 | Геометрический скрининг с водой 8CHO
- Цель: проверить влияние неподвижной кристаллографической воды на предскрининг.
- Скрипт/команда: `python enzyme-kinetics/scripts/screen_ksi_poses.py --wt enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --water-source enzyme-kinetics/data/research/p00947/8CHO-assembly.pdb --output enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water --seeds 101 202 303`.
- Входы: 8CHO WT SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`, 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`, исходная сборка 8CHO SHA-256 `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`.
- Параметры: RDKit 2026.09.1, жёсткое тело, искусственный геометрический скоринг, 48 атомов O кристаллографической воды по Bio.PDB, три seed, общая продуктовая матрица.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/poses-P00947-v3-water/summary.json` и три XYZ; журнал stdout `/tmp/opencode/ksi-poses-v3.log`.
- Результат: для 101/202/303 C4–O ближайшего H-пары 3.507/3.510/3.507 Å, H–O 2.439/2.839/2.439 Å, углы 164.19/119.38/164.19°. Скрининг расстояния прошёл для всех трёх, но эта величина была включена в целевую функцию.
- Отрицательный результат?: да — это не MM-энергия, не проверенный β-H и не подтверждённый WT–ES; ожидается совпадение части стартов, нужна количественная проверка различий.
- Диагностика: воды фиксированы, их каталитическая роль неизвестна; используются только тяжёлые атомы и эмпирические штрафы, белок и водная сеть не релаксировали.
- Влияние на следующие шаги: сравнить позы и подготовить физически параметризованные комплексы; для PMF/QM/MM результата недостаточно.
- Статус: partial.

## 2026-11-01 11:11 | Проверка кристаллографической воды
- Цель: проверить возможность переноса воды исходной сборки в нормализованный WT.
- Скрипт/команда: Bio.PDB: покомпонентное сравнение координат CA двух моделей исходной сборки и цепей A/B WT, подсчёт остатков HOH по моделям.
- Входы: `enzyme-kinetics/data/research/p00947/8CHO-assembly.pdb` SHA-256 `934adb04180d5dbc93ae88aafbb0d0967f307c11d9f51feb269110f2432164e3`; `enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb` SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`.
- Параметры: Bio.PDB, 125 CA на цепь, исходные координаты в Å.
- Выходы: вывод команды.
- Результат: 2 модели исходной биосборки; максимальное отклонение CA от A/B нормализованного WT 0.0 Å; по 24 воды на модель, прочитанные Bio.PDB.
- Отрицательный результат?: нет; каталитическая роль этих вод и сохранность всех waters с альтернативными положениями UNKNOWN.
- Диагностика: проверено совпадение системы координат; PDB-водороды отсутствуют, минимизация воды не выполнялась.
- Влияние на следующие шаги: можно отдельно повторить стерический скрининг с экспериментальными водами как фиксированными препятствиями, без объявления их каталитическими.
- Статус: done для проверки системы координат.

## 2026-11-01 11:09 | Повтор с проверкой воды
- Цель: исправить отчётность геометрического скрининга, не переписывая первый результат.
- Скрипт/команда: `python enzyme-kinetics/scripts/screen_ksi_poses.py --wt enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --output enzyme-kinetics/outputs/active-site-v3/poses-P00947-v2 --seeds 101 202 303`.
- Входы: WT SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`; 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`.
- Параметры: тот же скрининг и seeds, исправлен подсчёт входных атомов воды; RDKit 2026.09.1.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/poses-P00947-v2/summary.json` и три XYZ.
- Результат: рассчитано 0 атомов кристаллографической воды, флаг `experimental_water_included: false`; геометрические результаты прежние.
- Отрицательный результат?: да — воду, β-H и MM-энергию не установили; три старта не являются тремя различными позами.
- Диагностика: исходный безводный WT-файл не изменялся; никакая искусственная вода не добавлена.
- Влияние на следующие шаги: подготовить мотивированную водную сеть и физический скоринг; QM/MM по-прежнему не обоснован.
- Статус: partial.

## 2026-11-01 11:08 | Исправление к записи от 11:06: вода и различия поз
- Цель: проверить фактическое наличие воды и количественно сравнить координаты кандидатов.
- Скрипт/команда: Python stdlib, попарный RMSD тяжёлых атомов XYZ в неподвижной системе WT и подсчёт HETATM HOH/WAT в 8CHO-wt-dimer.pdb.
- Входы: `enzyme-kinetics/outputs/active-site-v3/poses-P00947/candidate-{101,202,303}.xyz`; WT SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`; SHA-256 XYZ см. вывод `sha256sum` в последующем воспроизводимом запуске.
- Параметры: 21 тяжёлый атом, без дополнительного выравнивания; анализ состава PDB.
- Выходы: вывод команды; исходные данные не изменялись.
- Результат: попарные RMSD 0.00061–0.00316 Å; в используемом WT PDB **0 молекул кристаллографической воды**.
- Отрицательный результат?: да — флаг `experimental_water_included: true` в исходном `summary.json` ошибочен и не должен использоваться; геометрический расчёт фактически не включал воду.
- Диагностика: флаг первоначально устанавливался без проверки PDB; теперь код подсчитывает воду явно. Первая запись и её артефакты сохранены как история.
- Влияние на следующие шаги: результаты не соответствуют требованию о водной сети и не являются независимым ансамблем; необходим новый запуск с исправленной отчётностью и последующая подготовка воды.
- Статус: done для выявления ошибки; water-aware исследование blocked.

## 2026-11-01 11:06 | Геометрический предскрининг трёх стартов
- Цель: проверить возможность трёх независимых стартов стероида у Asp38 в WT-димере.
- Скрипт/команда: `python enzyme-kinetics/scripts/screen_ksi_poses.py --wt enzyme-kinetics/outputs/experiments/active-site-v3/8CHO-wt-dimer.pdb --template enzyme-kinetics/data/research/p00947/3NHX-asymmetric.pdb --output enzyme-kinetics/outputs/active-site-v3/poses-P00947 --seeds 101 202 303`.
- Входы: WT 8CHO SHA-256 `2a7379e2297a97b6ebef3c81d49c2f5bcd636be5cac40080bbba83dc91b5ee13`; 3NHX SHA-256 `76cc10d1641569b45a48e7acfeab6455f23347c44103f0e2ee2f0eb8c082be5d`.
- Параметры: RDKit 2026.09.1, SciPy L-BFGS-B, жёсткое тело; разные seed, общий мутантный продуктовый шаблон, искусственная функция отталкивания и ограничения расстояния; не MM и не докинг.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/poses-P00947/` — три XYZ и `summary.json`.
- Результат: C4–O 3.184/3.181/3.185 Å; ближайший H–O 2.284/2.281/2.284 Å, C4–H–O 138.01/138.05/138.00°; все три старта сходятся практически к одной позе.
- Отрицательный результат?: да — трёх различных состояний или подтверждённой реакционной β-H-геометрии нет; энергия взаимодействия UNKNOWN.
- Диагностика: геометрия смещена искусственным ограничением и не доказывает связывание; протонирование, водная сеть и β-H не проверены; нет экспериментального WT–ES и независимых реплик.
- Влияние на следующие шаги: нужны альтернативные структуры/конформеры и физическая параметризация; эти кандидаты нельзя использовать для QM/MM как подтверждённый ES.
- Статус: partial — предварительная геометрия, а не ансамбль.

## 2026-11-01 10:50 | Проверка тестов химии
- Цель: проверить сопоставление структуры и устойчивость к неверным входам.
- Скрипт/команда: `python -m pytest -q`.
- Входы: `enzyme-kinetics/scripts/check_p00947_chemistry.py` SHA-256 `ded7bd569201d6b6b7e841839f659e28fcbc36dafb72609d7b4c7a706e99f6df`; `enzyme-kinetics/tests/test_p00947_chemistry.py` SHA-256 `9af3ed020196cf6adeb2571c718d0226679fe6983c3891c5c5d096f857e97d4b`.
- Параметры: Python 3.13, существующие site-packages, RDKit 2026.09.1.
- Выходы: вывод pytest.
- Результат: 55 passed за 5.52 s.
- Отрицательный результат?: нет.
- Диагностика: проверены правильный мотив, обратный порядок и отличие одного стереоцентра; тесты не подтверждают экспериментальный механизм.
- Влияние на следующие шаги: проверка химической записи воспроизводима.
- Статус: done.

## 2026-11-01 10:49 | Воспроизводимая проверка химии
- Цель: сопоставить локальную запись с зафиксированными эталонами PubChem и проверить реакционный граф.
- Скрипт/команда: `python enzyme-kinetics/scripts/check_p00947_chemistry.py --csv enzyme-kinetics/data/Baseline/kcat/kcat-random_trainval.csv --row 1526 --reference docs/experiments/active-site-v3/p00947-reference-molecules.json --output enzyme-kinetics/outputs/active-site-v3/p00947-chemistry-check.json`.
- Входы: CSV SHA-256 `99d6a214242dc709a5ab50475116386e1a204215565acf35888dce3621d8f44b`, сохранённый эталон SHA-256 `f1353cd685f1bf7de40310dcb2775cfd72816aa072cd86b80a99b83d992d2a29`.
- Параметры: RDKit 2026.09.1, стандартный InChIKey, canonical isomeric SMILES, граф исходного порядка атомов.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/p00947-chemistry-check.json`; код `enzyme-kinetics/scripts/check_p00947_chemistry.py`.
- Результат: оба сравнения SMILES и полных InChIKey совпадают, заряды 0; C4/C5/C6 в исходном порядке 6/7/8, кратности C4–C5 и C5–C6 изменяются 1/2 → 2/1.
- Отрицательный результат?: нет; β-H и точная стереохимия партии assay ещё не установлены.
- Диагностика: проверены 21 атом, рёбра и элементы в исходных графах S/P; индексы применимы только к этой записи; равенство связности не доказывает механизм.
- Влияние на следующие шаги: можно использовать подтверждённую химическую идентичность PubChem как основу подготовки поз, не считая её подтверждением WT–ES.
- Статус: done для сравнения записей, partial для стереоспецифического переноса.

## 2026-11-01 10:47 | Проверка атомных индексов
- Цель: проверить реакционную нумерацию для исходного SMILES.
- Скрипт/команда: RDKit `MolFromSmiles`, вывод индексов, соседей и кратности связей локального субстрата.
- Входы: локальный SMILES из `reaction-card-P00947.md`, данные исходной CSV SHA-256 `99d6a214242dc709a5ab50475116386e1a204215565acf35888dce3621d8f44b`.
- Параметры: RDKit 2026.09.1, нумерация атомов от 0 и порядок исходного SMILES.
- Выходы: вывод команды; воспроизводимый контроль запланирован отдельно.
- Результат: в этой записи мотив C3–C4–C5=C6 соответствует индексам 4–6–7–8; другая индексация 13/15 в старом кластере относится к другому порядку атомов и не переносится автоматически.
- Отрицательный результат?: нет; стереохимическое назначение конкретного β-H остаётся UNKNOWN.
- Диагностика: C3 идентифицирован по соседнему O карбонила, C5 — разветвлённый атом двойной связи внутри шестичленного кольца; требуется проверка сопоставления с продуктом.
- Влияние на следующие шаги: при подготовке структуры нельзя переносить старые индексы без атомного сопоставления.
- Статус: partial.

## 2026-11-01 10:46 | Первичное сравнение стереохимии
- Цель: сопоставить локальные субстрат и продукт с независимыми химическими записями.
- Скрипт/команда: RDKit `MolToSmiles`, `CalcMolFormula`, `MolToInchiKey`, `FindMolChiralCenters` для локальных строк и PubChem SMILES; первичная интерактивная проверка.
- Входы: `enzyme-kinetics/data/Baseline/kcat/kcat-random_trainval.csv` SHA-256 `99d6a214242dc709a5ab50475116386e1a204215565acf35888dce3621d8f44b`; PubChem PUG REST CID 160531 и 6128, URL `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/160531,6128/property/IsomericSMILES,InChI,InChIKey,MolecularFormula/JSON`; содержимое веб-ответа в этом запуске не сохранено, SHA-256 UNKNOWN.
- Параметры: RDKit 2026.09.1, изомерные SMILES и стандартные InChIKey.
- Выходы: пока только вывод команды; воспроизводимая фиксация будет отдельно.
- Результат: локальный субстрат и CID 160531 — `SQGZFRITSMYKRH-QAGGRKNESA-N`; локальный продукт и CID 6128 — `AEMFNILZOJDQLW-QAGGRKNESA-N`; все четыре C19H26O2, пять заданных тетраэдрических центров.
- Отрицательный результат?: нет; полное соответствие конкретной партии реактива статье и позиционная нумерация реакционных атомов этим не доказываются.
- Диагностика: отдельно проверены канонические изомерные SMILES и стандартные InChIKey; сравнение неподтверждённой индексной нумерации атомов не использовалось.
- Влияние на следующие шаги: необходимы фиксированный снимок эталонов, сопоставление атомов и тесты; WT–ES остаётся неподтверждённым.
- Статус: partial.

## 2026-11-01 10:10 | Проверка страницы методов
- Цель: уточнить концентрации и метод assay по первичному источнику.
- Скрипт/команда: `pdftotext -f 2 -l 2 -layout <PDF> -` — дополнительная проверка после преобразования PDF в Markdown.
- Входы: `articles/Origin of the Different pH Activity Profile in Two Homologous Ketosteroid Isomerases*.pdf`, SHA-256 `bf405ad64046bacc66070ddae8793b4bae8254be2408ec125046da45224d70a0`.
- Параметры: текстовый слой, локальное извлечение; символ микро отображён заменяющим знаком, единица проверяется в контексте.
- Выходы: уточнение `reaction-card-P00947.md`; файл pdftotext не сохранялся.
- Результат: 5-AND 12, 35, 58, 82, 116 мкМ, метанол 3.3 об.%; начальная скорость за 1–2 мин, превращение <10%, Lineweaver–Burk.
- Отрицательный результат?: нет; природа ошибок ± и полная стереохимия всё ещё UNKNOWN.
- Диагностика: сохранены раздельными условия основного assay и pH-серии; значения взяты из раздела методов, не из подписи графика.
- Влияние на следующие шаги: уточнить карточку; не считать условия автоматически тождественными исходным CSV.
- Статус: partial — биохимическое тождество локальной записи всё ещё проверяется.

## 2026-11-01 10:09 | Локальная статья и MarkItDown
- Цель: получить машиночитаемый текст первичной статьи и проверить assay.
- Скрипт/команда: `uv tool install 'markitdown[pdf]'`; `markitdown 'articles/Origin of the Different pH Activity Profile in Two Homologous Ketosteroid Isomerases*.pdf' -o articles/Yun-2003-pH-profile.md`; дополнительная сверка таблицы: `pdftotext -f 3 -l 3 -layout <PDF> -`.
- Входы: `articles/Origin of the Different pH Activity Profile in Two Homologous Ketosteroid Isomerases*.pdf`, SHA-256 `bf405ad64046bacc66070ddae8793b4bae8254be2408ec125046da45224d70a0`.
- Параметры: markitdown 0.1.9, PDF extra, изолированная пользовательская установка uv; плагины и облачные сервисы не применялись.
- Выходы: `articles/Yun-2003-pH-profile.md`, SHA-256 `27f5a68f8fb7afd55fc184d5466e0feab60d272b185eac70083269f0d44ee263`; обновлённая карточка.
- Результат: таблица I, с. 28231: TI-WT kcat 31683 ± 1200 с⁻¹, Km 152.8 ± 8.2 мкМ; условия assay включают 3.3 об.% метанола, 5 концентраций 5-AND.
- Отрицательный результат?: да — Markdown искажает двухколоночный порядок и отдельные символы; значение ± не установлено по проверенным фрагментам, стереохимия локального SMILES независимо не доказана.
- Диагностика: заголовки единиц и выравнивание таблицы дополнительно проверены локальным `pdftotext -layout`, без выгрузки в сеть; не смешивать таблицу TI-WT с PI-WT и pKa-серией с другим субстратом.
- Влияние на следующие шаги: карточку можно уточнить; для WT–ES всё ещё необходимо сверить структуру субстрата и геометрию.
- Статус: partial — текст и таблица доступны; полное стереохимическое соответствие и природа ошибок остаются неизвестными.

## 2026-11-01 09:56 | Отбор поз и граница QM/MM
- Цель: оценить готовность известных входов для WT–ES и условного QM/MM-пилота.
- Скрипт/команда: чтение сохранённых `H2-conditions-factorial.json`, `H5-mechanism-status.json`, `H5-alternative-product-template.json`; новых PDB2PQR, докинга, минимизации или QM/MM не запускалось.
- Входы: H2 SHA-256 `19759155571f72272b9acaef7a3b83b18d8903520d16dafff4c05f9b7e5b7ed3`; H5 SHA-256 `1299164eacde4637764489af89134715f418f511757a23fdf6d80bf2999bc98d`; H5 alternative SHA-256 `9b2a76aad12228cf4b2753c66dbcd29e81573ba6ad306e8ae771e09c2454bc42`.
- Параметры: старые apo PROPKA и геометрические пробы; без новых seed.
- Выходы: [poses-P00947.md](poses-P00947.md), [mechanisms-P00947.md](mechanisms-P00947.md).
- Результат: прежняя 1OHS-поза отклонена, альтернативная 3NHX не валидирована как WT–ES; 3 новых кандидата не получены.
- Отрицательный результат?: да — подтверждённого комплекса нет; условие запуска QM/MM не выполнено.
- Диагностика: локальная реакция и первичный assay полностью по стереохимии не сопоставлены; продуктовый шаблон мутантный, атом H не верифицирован.
- Влияние на следующие шаги: требуется первичная статья/подтверждённая структура субстрата и новая подготовка поз; невозможно заключить, что все позы невозможны.
- Статус: blocked — неподтверждённые входы; новые гипотезы являются протоколом, а не расчётными результатами.

## 2026-11-01 09:54 | Полный набор тестов
- Цель: проверить проект после добавления диагностики.
- Скрипт/команда: `python -m pytest -q`.
- Входы: тесты `enzyme-kinetics/tests/`; изменённый `test_md_rmsd_consistency.py` SHA-256 `7cdc5cdee82c5672f5da70ffece4c33d834fe8d8ffb72a80cbc25da637779e7d`, скрипт диагностики SHA-256 `56ec7f24d553cefd327509f1cdf29c888db8b9b8c5f9cc036080edb777d0262f`.
- Параметры: Python 3.13, site-packages существующего окружения; MDTraj как необязательная зависимость.
- Выходы: вывод pytest.
- Результат: 52 passed за 5.91 s.
- Отрицательный результат?: нет.
- Диагностика: полный локальный pytest; научную корректность поз и реакций тесты не доказывают.
- Влияние на следующие шаги: программный набор прошёл, подтверждение ES и assay остаётся открытым.
- Статус: done.

## 2026-11-01 09:52 | Тесты диагностики и аудита
- Цель: проверить восстановление периодического образа и регрессии аудита.
- Скрипт/команда: `python -m pytest -q enzyme-kinetics/tests/test_md_rmsd_consistency.py enzyme-kinetics/tests/test_reaction_registry_audit.py`.
- Входы: `enzyme-kinetics/scripts/diagnose_dimer_rmsd.py` SHA-256 `56ec7f24d553cefd327509f1cdf29c888db8b9b8c5f9cc036080edb777d0262f`; `enzyme-kinetics/tests/test_md_rmsd_consistency.py` SHA-256 `0a59b8ac4f70200546309fdb0ab944786dac88461072dc627e389335db19f8a5`; `enzyme-kinetics/tests/test_reaction_registry_audit.py` SHA-256 `c64814e2a37ffc4c00d2e2667826aab6d0c05b6df5e10c92b4d18aaafb1f0b89`.
- Параметры: Python 3.13, установленный отдельно; site-packages существующего окружения.
- Выходы: вывод pytest; новые файлы данных отсутствуют.
- Результат: 6 passed.
- Отрицательный результат?: нет.
- Диагностика: синтетический скачок на размер ячейки устраняется, нулевая ячейка отвергается.
- Влияние на следующие шаги: тесты подтверждают узкий диагностический сценарий, но не весь конвейер.
- Статус: done.

## 2026-11-01 09:51 | Replica 101: периодические образы
- Цель: установить причину скачка RMSD.
- Скрипт/команда: `enzyme-kinetics/scripts/diagnose_dimer_rmsd.py --trajectory enzyme-kinetics/outputs/experiments/active-site-v3/replica-101/trajectory.dcd --topology enzyme-kinetics/outputs/experiments/active-site-v3/replica-101/topology.pdb --output enzyme-kinetics/outputs/active-site-v3/rmsd-replica-101-diagnostic.json`.
- Входы: DCD SHA-256 `2a8b5fdf51804727b726803e48c58edea2e232325477e9a42e5277d1bad5d618`; topology SHA-256 `db4d3c954e51ea7067811b4774aea5fa52102399407998749c911ba1a641d10b`.
- Параметры: MDTraj 1.11.2, 250 CA, 50 кадров, ортогональный box 7.903 нм; одна и та же референсная структура кадр 0; перенос цепи B по minimum image, без изменения траектории.
- Выходы: `enzyme-kinetics/outputs/active-site-v3/rmsd-replica-101-diagnostic.json`.
- Результат: исходный RMSD 0.05359/1.23276 нм по половинам воспроизведён; в 10 кадрах цепь B в соседнем периодическом образе на 7.903 нм. После реконструкции RMSD 0.05359/0.07075 нм; центры цепей 2.414–2.438 нм.
- Отрицательный результат?: да — прежний анализ не реконструировал периодический димер, поэтому его RMSD/RMSF/DCCM для replica 101 нельзя использовать без перерасчёта.
- Диагностика: один набор CA и кадр 0 для обеих половин, раздельные RMSD цепей не показывают аналогичного скачка; физическая стабильность не доказана.
- Влияние на следующие шаги: требуется корректный PBC-перерасчёт всей серии и её статистики; диагностика не изменяет исходные файлы.
- Статус: done для причины скачка; реконструкция предполагает целостность каждой цепи и не является оценкой сходимости.

## 2026-11-01 09:47 | Первичный поиск assay
- Цель: установить публикацию и условия для P00947.
- Скрипт/команда: поиск по DOI 10.1074/jbc.M302166200, чтение PubMed 12734184; попытка открыть полную статью ScienceDirect.
- Входы: `enzyme-kinetics/data/Baseline/kcat/kcat-random_trainval.csv` SHA-256 `99d6a214242dc709a5ab50475116386e1a204215565acf35888dce3621d8f44b`; `enzyme-kinetics/data/Baseline/km/km-random_trainval.csv` SHA-256 `d94608500a0fa049e0082f1eeb4b12603f56ee5b4a540f1772666557d8a524ec`; веб-страницы без локальной копии, SHA-256 UNKNOWN.
- Параметры: полнотекстовый доступ через webfetch; индексированный поисковый фрагмент отдельно от открытого PubMed.
- Выходы: [reaction-card-P00947.md](reaction-card-P00947.md).
- Результат: идентифицированы Yun et al. 2003, PMID 12734184; индекс содержит 5-AND, 34 mM фосфата калия, 2.5 mM EDTA, pH 7, 25 °C и детекцию 248 nm.
- Отрицательный результат?: да — полный текст ScienceDirect и Academia вернули HTTP 403; таблица и её заголовки непосредственно не прочитаны.
- Диагностика: PubMed подтверждает идентификаторы и организм, но не табличные единицы; индексный фрагмент не заменяет исходную статью.
- Влияние на следующие шаги: отсутствующие единицы, ошибки и точное соответствие локальных строк assay требуют первичной проверки.
- Статус: partial — требуется доступ к первичной таблице и проверка идентичности молекул.
