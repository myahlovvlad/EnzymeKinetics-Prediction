"""Export portable physical-model benchmark evidence without copying datasets/weights."""

import argparse
import json
from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import digest_file, write_json


def export(runs, output):
    output.mkdir(parents=True, exist_ok=True)
    rows, summaries = [], []
    fig, axes = plt.subplots(
        len(runs), 1, figsize=(10, max(4, 3 * len(runs))), squeeze=False
    )
    keep = [
        "target",
        "target_unit",
        "selected_model",
        "seed",
        "rows",
        "split_scheme",
        "group_by",
        "ec_level",
        "test_group_overlap_rows",
        "canonical_pair_overlap_test_rows",
        "feature_version",
        "feature_blocks",
        "physical_coverage",
        "coverage",
        "versions",
        "embeddings",
        "trainval_sha256",
        "test_sha256",
        "limitations",
        "elapsed_seconds",
        "code_sha256_at_start",
    ]
    for index, run in enumerate(runs):
        meta = json.loads((run / "metadata.json").read_text())
        metrics = json.loads((run / "metrics.json").read_text())
        destination = output / run.name
        destination.mkdir(exist_ok=True)
        shutil.copyfile(run / "report.md", destination / "README.md")
        shutil.copyfile(run / "metrics.json", destination / "metrics.json")
        summary = {k: meta[k] for k in keep}
        summary["artifact_checksums"] = {
            f: digest_file(run / f)
            for f in [
                "metadata.json",
                "metrics.json",
                "predictions.csv",
                "preprocessing.joblib",
                "split.csv",
            ]
        }
        write_json(destination / "summary.json", summary)
        selected = meta["selected_model"]
        for name, result in metrics.items():
            rows.append(
                {
                    "run": run.name,
                    "model": name,
                    "selected_by_validation": name == selected,
                    "validation_mae": result["validation"]["mae"],
                    **{
                        "test_pair_disjoint_" + k: v
                        for k, v in result["test_pair_disjoint"].items()
                    },
                }
            )
        summaries.append((run.name, meta, metrics[selected]["test_pair_disjoint"]))
        names = list(metrics)
        x = np.arange(len(names))
        axis = axes[index, 0]
        axis.bar(
            x - 0.18,
            [metrics[n]["validation"]["mae"] for n in names],
            0.36,
            label="Validation",
        )
        axis.bar(
            x + 0.18,
            [metrics[n]["test_pair_disjoint"]["mae"] for n in names],
            0.36,
            label="Test, pair disjoint",
        )
        axis.set_xticks(x, [n + ("*" if n == selected else "") for n in names])
        axis.set_ylabel("MAE, log10")
        axis.set_title(run.name + " (* selected on validation)")
        axis.legend(loc="upper right")
        axis.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(output / "comparison.png", dpi=160)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(output / "comparison.csv", index=False)
    lines = [
        "# Физические признаки и хемометрика: эксперименты v2",
        "",
        "Дата: 1 ноября 2026. Один seed: 42. Модель выбирается по validation.",
        "",
        "| Run | Модель | Train / val / test | MAE test без совпадающих пар | R² | Пересечение test-групп |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for name, meta, result in summaries:
        counts = meta["rows"]
        lines.append(
            f"| [{name}]({name}/README.md) | {meta['selected_model']} | {counts['train']} / {counts['validation']} / {counts['test']} | {result['mae']:.4f} | {result['r2']:.4f} | {meta['test_group_overlap_rows']} |"
        )
    lines += [
        "",
        "![Сравнение моделей](comparison.png)",
        "",
        "У provided-test сохраняется исходный test, а validation разнесена по exact sequence. Поэтому пересечение test-групп здесь ожидаемо и отражает ограничение исходного split.",
        "",
        "EC-OOD и cluster40 используют новые grouped-three-way части. Нулевое пересечение EC/cluster ID не доказывает аудит происхождения кластеров и отсутствие иных видов dataset leakage. Числа разных split нельзя напрямую ранжировать как улучшение модели.",
        "",
        "По сравнению с v1 kcat MLP ухудшилась 0.9609 → 0.9883; Km MLP немного улучшилась 0.6853 → 0.6797. Один seed и низкое покрытие физических признаков не подтверждают устойчивую пользу. Старый test уже просматривался; это исследовательские сравнения, нужен новый внешний holdout.",
        "",
        "3 841 pocket-запись; 34 ошибки extraction. QM: 4 успеха из 32 выбранных SMILES, 28 отказов преимущественно из-за разъединённых реагентов. Покрытие по train/validation/test записано в summary.json. Электростатические и MD-таблицы пилота не включены в массовое обучение: неизвестная ионная сила assay не заменяется предполагаемой.",
        "",
        "[Реальный пилот физических backend](physical-pilot/manifest.json) и его артефакты не являются доказательством свободноэнергетического барьера. Production MD — 0.2 ps, недостаточная термализация. ESM-2 рассчитан на двух последовательностях, не на всём датасете. [Статусы и следующий цикл](../../IMPLEMENTATION_STATUS.md).",
        "",
        "[Полная таблица](comparison.csv), [проверка reload](reload-verification.json). В подпапках — все метрики и SHA-256 локальных артефактов.",
    ]
    (output / "README.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.runs, args.output)


if __name__ == "__main__":
    main()
