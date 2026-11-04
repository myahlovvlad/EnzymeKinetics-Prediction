# Changelog

## 0.3.0 — 2026-11-01

- WT-димер 8CHO, каталитические роли, контактные остатки, геометрия chi1/chi2.
- Sequence pI/Q(pH)/pH−pI в новых признаках; совместимость старых checkpoint.
- Зарегистрированные гипотезы, журнал, дневник, факторные APBS-эксперименты,
  проверка сеток и focused-grid refinement; векторизованный OpenDX reader.
- Три независимые MD-реплики, RMSF/DCCM/ротамеры, exact replicate-level tests.
- Условные торсионные потенциалы и неаддитивность Tyr14/Asp99.
- Явная ES/TS1/EI/TS2/EP-схема, проверка Hessian; неподходящая поза отклонена.
- Измеренный бюджет CPU 3×8, проверка GPU; MD checkpoint continuation.
- Переобученная kcat MLP с титрационными признаками; научные ограничения
  и отрицательные результаты сохранены в docs/experiments.

## 0.2.0 — 2026-11-01

- Реальные P2Rank, GFN2-xTB, PDB2PQR, OpenMM; реальный APBS LPBE с OpenDX parser
  и установкой официальной сборки по закреплённой SHA-256.
- Артефакты с единицами, provenance, checksum, конфигурацией и статусами.
- Физические таблицы в train/inference, Consensus MB-PLS, single-y OPLS,
  EC-local проекции с глобальным fallback; EC/cluster grouping.
- ESM-2 exporter с immutable revision и явной политикой длинных последовательностей.
- Обратимая mass-action реакция, conservation/Haldane checks, TST/tunneling.
- Анализ all-atom trajectory distances и проекции приближённого поля.
- Реальные kcat/Km и EC/cluster эксперименты; native tests и CI matrix.
- Совместимость существующих v1 weights/preprocessing сохранена.

Научные ограничения и неподтверждённые части — docs/IMPLEMENTATION_STATUS.md.

## 0.1.0 — 2026-11-01

- Единый корневой uv-проект, перенос уникальных материалов и аудит дублей.
- Воспроизводимое обучение kcat/Km, validation-only selection и inference.
- Базовые линейные, латентные и нейросетевые модели; сохранённые v1-отчёты.
