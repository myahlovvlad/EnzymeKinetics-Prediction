"""Generate actual ESM-2 residue-mean embeddings with a recorded model revision."""

import argparse
from pathlib import Path

import pandas as pd
import torch

from enzymekinetics_prediction.physics.artifacts import (
    digest_file,
    package_version,
    write_json,
)


def embed_sequences(
    sequences,
    output,
    model_id="facebook/esm2_t6_8M_UR50D",
    revision=None,
    batch_size=4,
    threads=2,
):
    if not model_id.startswith("facebook/esm2_"):
        raise ValueError("This provider implements official facebook ESM-2 models only")
    if batch_size < 1 or threads < 1:
        raise ValueError("Positive batch and thread counts required")
    try:
        from huggingface_hub import model_info
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:
        raise RuntimeError("Install the embeddings extra for ESM-2") from exc
    sequences = sorted(set(str(s).strip() for s in sequences))
    if not sequences or any(
        not s or any(a not in "ACDEFGHIKLMNPQRSTVWYXBZUO" for a in s) for s in sequences
    ):
        raise ValueError("Invalid or empty protein sequence")
    if revision is None:
        revision = model_info(model_id).sha
    if len(revision) != 40 or any(
        c not in "0123456789abcdef" for c in revision.lower()
    ):
        raise ValueError("Use an immutable 40-character model revision")
    output = Path(output)
    if output.exists():
        raise ValueError("Embedding output already exists")
    tokenizer = AutoTokenizer.from_pretrained(
        model_id, revision=revision, trust_remote_code=False
    )
    model = AutoModel.from_pretrained(
        model_id, revision=revision, trust_remote_code=False, add_pooling_layer=False
    )
    model.eval().cpu()
    max_residues = int(model.config.max_position_embeddings) - 2
    if any(len(s) > max_residues for s in sequences):
        raise ValueError(
            f"Long sequences require an explicit chunking protocol; maximum here is {max_residues}"
        )
    torch.set_num_threads(threads)
    vectors = {}
    with torch.no_grad():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]
            tokens = tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                truncation=False,
                return_special_tokens_mask=True,
            )
            special = tokens.pop("special_tokens_mask").bool()
            output_hidden = model(**tokens).last_hidden_state
            mask = tokens["attention_mask"].bool() & ~special
            pooled = (output_hidden * mask.unsqueeze(-1)).sum(dim=1) / mask.sum(
                dim=1, keepdim=True
            )
            for sequence, vector in zip(batch, pooled, strict=True):
                vectors[sequence] = vector.detach().cpu().float()
    output.parent.mkdir(parents=True, exist_ok=True)
    torch.save(vectors, output)
    metadata = {
        "provider": "ESM-2",
        "model_id": model_id,
        "revision": revision,
        "pooling": "residue mean excluding special/padding tokens",
        "sequences": len(vectors),
        "dimension": len(next(iter(vectors.values()))),
        "sha256": digest_file(output),
        "long_sequence_policy": "reject, no silent truncation",
        "sequence_keys": "exact sequence strings",
        "versions": {
            p: package_version(p) for p in ["torch", "transformers", "huggingface-hub"]
        },
    }
    write_json(output.with_suffix(output.suffix + ".meta.json"), metadata)
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model-id", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--revision")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    frame = pd.read_csv(args.input)
    metadata = embed_sequences(
        frame.sequence.dropna(),
        args.output,
        args.model_id,
        args.revision,
        args.batch_size,
        args.threads,
    )
    print(metadata)


if __name__ == "__main__":
    main()
