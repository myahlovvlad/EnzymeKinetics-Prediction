# Рабочая кодовая база

В этом каталоге находятся пакет, скрипты, тесты, исследовательские документы
и локальные данные. Зависимости и Git управляются из родительского корня.
Отдельных `pyproject.toml`, `uv.lock` и `.venv` здесь нет.

Из этого каталога `uv` автоматически находит корневой проект:

```sh
uv run --locked enzyme-train --target kcat
uv run --locked pytest -q
uv run --locked python scripts/diagnostics/check_environment.py --profile chemistry --profile quantum --profile dashboard --profile ml
```

- [Главный README](../README.md)
- [Обучение и inference](../docs/TRAINING.md)
- [Результаты kcat](../docs/benchmarks/kcat/README.md)
- [Результаты Km](../docs/benchmarks/km/README.md)
- [Исследовательская постановка](docx/README.md)

`apps/catalysis.py` — иллюстративная анимация; она не использует обученные веса
и не выполняет молекулярную динамику. Предсказания выполняет `enzyme-predict`.

Физический контур версии 0.2: [инструкция](../docs/PHYSICS.md).
Соответствие исходной идее и критерии следующего этапа:
[матрица требований](../docs/IMPLEMENTATION_STATUS.md).
