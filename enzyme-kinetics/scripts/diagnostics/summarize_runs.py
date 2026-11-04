"""Export portable benchmark reports and figures from completed local runs."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from enzymekinetics_prediction.learning.features import enzyme_substrate_pairs
from enzymekinetics_prediction.learning.training import metrics
from enzymekinetics_prediction.paths import DATA_DIR


def export_run(run, output):
    metadata = json.loads((run / "metadata.json").read_text())
    target = metadata["target"]
    scores = json.loads((run / "metrics.json").read_text())
    predictions = pd.read_csv(run / "predictions.csv")
    trainval = pd.read_csv(DATA_DIR / f"Baseline/{target}/{target}-random_trainval.csv")
    test = pd.read_csv(DATA_DIR / f"Baseline/{target}/{target}-random_test.csv")
    if len(test) != len(predictions):
        raise ValueError(
            "Prediction rows do not match dataset; audit valid source records first"
        )
    pairs = set(enzyme_substrate_pairs(trainval))
    clean = np.asarray([pair not in pairs for pair in enzyme_substrate_pairs(test)])
    predictions["canonical_pair_disjoint_vs_trainval"] = clean
    for name in scores:
        scores[name]["test_pair_disjoint"] = metrics(
            predictions.log10_value[clean], predictions[name][clean]
        )
    (run / "metrics.json").write_text(json.dumps(scores, indent=2) + "\n")
    predictions.to_csv(run / "predictions.csv", index=False)
    destination = output / target
    destination.mkdir(parents=True, exist_ok=True)
    summary = {
        k: metadata[k]
        for k in [
            "target",
            "target_column",
            "target_scale",
            "target_unit",
            "seed",
            "threads",
            "selected_model",
            "selection_metric",
            "rows",
            "coverage",
            "versions",
            "feature_version",
            "trainval_sha256",
            "test_sha256",
            "elapsed_seconds",
            "limitations",
        ]
    }
    summary["embedding_info"] = {
        k: v for k, v in metadata["embeddings"].items() if k != "path"
    }
    summary["canonical_pair_overlap_test_rows"] = int((~clean).sum())
    summary["neural_epochs"] = {
        k: {a: metadata["models"][k][a] for a in ["best_epoch", "epochs_run"]}
        for k in ["mlp", "radial"]
        if k in metadata["models"]
    }
    summary["code_provenance"] = (
        "v1 run completed before source commits; data hashes, split, weights and versions preserved"
    )
    (destination / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n"
    )
    (destination / "metrics.json").write_text(json.dumps(scores, indent=2) + "\n")
    lines = [
        f"# {target}: завершённое обучение",
        "",
        f"Выбор по validation: **{metadata['selected_model']}**. Seed: {metadata['seed']}.",
        "",
        "| Модель | Validation MAE | Test MAE | Test без совпадающих пар: MAE | R² на очищенном test | MAE для новых точных последовательностей |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for name, score in scores.items():
        lines.append(
            f"| {name} | {score['validation']['mae']:.4f} | {score['test']['mae']:.4f} | {score['test_pair_disjoint']['mae']:.4f} | {score['test_pair_disjoint']['r2']:.4f} | {score['test_unseen_exact_sequence']['mae']:.4f} |"
        )
    lines += [
        "",
        f"Размеры: {metadata['rows']}. Исключено совпадающих канонических пар: {int((~clean).sum())}.",
        "Все MAE в log10. Очистка test выполнена после обучения без переобучения и без смены выбранной модели.",
        "Новые точные последовательности не означают независимость от гомологов. Это один seed, а не оценка статистической значимости.",
        "Сводки геометрии доступны примерно для половины записей; доли пропусков приведены в summary.json.",
        "",
        "![Validation и test](comparison.png)",
        "",
        "![Предсказания выбранной модели](predictions.png)",
    ]
    report = "\n".join(lines) + "\n"
    (destination / "README.md").write_text(report)
    # Keep the local run report consistent with the audited, portable table.
    (run / "report.md").write_text(
        report.replace(
            "(comparison.png)",
            "(../../../../docs/benchmarks/" + target + "/comparison.png)",
        ).replace(
            "(predictions.png)",
            "(../../../../docs/benchmarks/" + target + "/predictions.png)",
        )
    )
    names = list(scores)
    fig, ax = plt.subplots(figsize=(8, 4))
    positions = np.arange(len(names))
    ax.bar(
        positions - 0.18,
        [scores[n]["validation"]["mae"] for n in names],
        width=0.36,
        label="Grouped validation",
    )
    ax.bar(
        positions + 0.18,
        [scores[n]["test_pair_disjoint"]["mae"] for n in names],
        width=0.36,
        label="Pair-disjoint test",
    )
    ax.set_xticks(positions, names)
    ax.set_ylim(0, 1.35 * max(scores[n]["validation"]["mae"] for n in names))
    ax.set_ylabel("MAE (log10 target)")
    ax.set_title(f"{target}: one-seed local benchmark")
    ax.legend()
    fig.tight_layout()
    fig.savefig(destination / "comparison.png", dpi=160)
    plt.close(fig)
    selected = metadata["selected_model"]
    y, predicted = (
        predictions.log10_value.to_numpy()[clean],
        predictions[selected].to_numpy()[clean],
    )
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.hexbin(y, predicted, gridsize=35, mincnt=1, cmap="viridis")
    limits = [
        float(min(y.min(), predicted.min())),
        float(max(y.max(), predicted.max())),
    ]
    ax.plot(limits, limits, color="red", linestyle="--", linewidth=1)
    ax.set_xlabel("Observed (log10 dataset units)")
    ax.set_ylabel("Predicted (log10 dataset units)")
    ax.set_title(f"{target}: {selected}, pair-disjoint test")
    fig.tight_layout()
    fig.savefig(destination / "predictions.png", dpi=160)
    plt.close(fig)
    print(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for run in args.runs:
        export_run(run, args.output)


if __name__ == "__main__":
    main()
