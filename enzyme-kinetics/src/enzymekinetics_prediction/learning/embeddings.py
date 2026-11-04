"""Generate actual ESM-2 residue-mean embeddings with a recorded model revision."""

import argparse
import json
from pathlib import Path

import pandas as pd
import torch

from enzymekinetics_prediction.physics.artifacts import (
    digest_file,
    package_version,
    write_json,
    content_key,
)

ENCODING_ALGORITHM = "esm2-residue-mean-excluding-special-padding-v1"
SOURCE_FILES = (
    "model.safetensors",
    "config.json",
    "vocab.txt",
    "tokenizer_config.json",
    "special_tokens_map.json",
)


def cached_sources(model_id, revision):
    """Inspect local source files only; never resolve or download missing files."""
    from huggingface_hub import try_to_load_from_cache

    return {
        name: digest_file(path)
        for name in SOURCE_FILES
        if isinstance(
            path := try_to_load_from_cache(model_id, name, revision=revision), str
        )
        and Path(path).is_file()
    }


def embedding_identity(sequences, model_id, revision, sources):
    """Content identity for strict H1 caches; locations are not identity."""
    if set(sources) != set(SOURCE_FILES):
        raise ValueError(
            "Strict embedding cache requires all model/tokenizer source hashes"
        )
    sequences = sorted(set(str(s).strip() for s in sequences))
    if (
        not sequences
        or not isinstance(revision, str)
        or len(revision) != 40
        or any(c not in "0123456789abcdef" for c in revision)
    ):
        raise ValueError(
            "Embedding identity requires sequences and an immutable revision"
        )
    if any(
        not isinstance(v, str)
        or len(v) != 64
        or any(c not in "0123456789abcdef" for c in v)
        for v in sources.values()
    ):
        raise ValueError("Embedding identity requires SHA-256 source digests")
    return {
        "schema_version": 1,
        "algorithm": ENCODING_ALGORITHM,
        "algorithm_code_sha256": digest_file(__file__),
        "model_id": model_id,
        "revision": revision,
        "source_sha256": dict(sources),
        "sequence_sha256": content_key(sequences),
        "sequence_count": len(sequences),
        "sequence_normalization": "strip, sorted unique; reject invalid; no truncation",
        "pooling": "residue mean excluding special/padding tokens",
        "output_dtype": "float32",
        "frozen": True,
        "versions": {
            p: package_version(p) for p in ["torch", "transformers", "huggingface-hub"]
        },
    }


def validate_embedding_cache(path, expected_identity, dimension):
    """Reject stale provenance or altered payloads before strict H1 training."""
    path = Path(path)
    metadata = json.loads(path.with_suffix(path.suffix + ".meta.json").read_text())
    if metadata.get("identity") != expected_identity or metadata.get(
        "cache_key"
    ) != content_key(expected_identity):
        raise ValueError(
            "Embedding cache identity mismatch or legacy unversioned sidecar"
        )
    if metadata.get("sha256") != digest_file(path):
        raise ValueError("Embedding cache payload checksum mismatch")
    vectors = torch.load(path, map_location="cpu", weights_only=True)
    if (
        not isinstance(vectors, dict)
        or not vectors
        or any(not isinstance(k, str) for k in vectors)
    ):
        raise ValueError("Embedding cache requires exact sequence keys")
    if (
        content_key(sorted(vectors)) != expected_identity["sequence_sha256"]
        or len(vectors) != expected_identity["sequence_count"]
    ):
        raise ValueError("Embedding cache sequence membership mismatch")
    if metadata.get("dimension") != dimension or any(
        not isinstance(v, torch.Tensor)
        or v.shape != (dimension,)
        or v.dtype != torch.float32
        or not torch.isfinite(v).all()
        for v in vectors.values()
    ):
        raise ValueError("Embedding cache dimension/dtype/finite-value mismatch")
    return metadata["cache_key"]


def embed_sequences(
    sequences,
    output,
    model_id="facebook/esm2_t6_8M_UR50D",
    revision=None,
    batch_size=4,
    threads=2,
    local_files_only=False,
):
    if not model_id.startswith("facebook/esm2_"):
        raise ValueError("This provider implements official facebook ESM-2 models only")
    if batch_size < 1 or threads < 1:
        raise ValueError("Positive batch and thread counts required")
    if local_files_only and revision is None:
        raise ValueError("Offline encoding requires an explicit immutable revision")
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
    if output.exists() or output.with_suffix(output.suffix + ".meta.json").exists():
        raise ValueError("Embedding output already exists")
    tokenizer = AutoTokenizer.from_pretrained(
        model_id,
        revision=revision,
        trust_remote_code=False,
        local_files_only=local_files_only,
    )
    model = AutoModel.from_pretrained(
        model_id,
        revision=revision,
        trust_remote_code=False,
        add_pooling_layer=False,
        local_files_only=local_files_only,
        use_safetensors=True if local_files_only else None,
    )
    model.eval().cpu()
    model.requires_grad_(False)
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
    sources = cached_sources(model_id, revision)
    identity = (
        embedding_identity(sequences, model_id, revision, sources)
        if set(sources) == set(SOURCE_FILES)
        else None
    )
    if local_files_only and identity is None:
        raise ValueError("Offline H1 cache has incomplete tokenizer/model provenance")
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
        "local_files_only": local_files_only,
        "frozen": True,
        "cached_source_sha256": sources,
        "identity": identity,
        "cache_key": content_key(identity) if identity is not None else None,
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
    parser.add_argument("--local-files-only", action="store_true")
    args = parser.parse_args()
    frame = pd.read_csv(args.input)
    metadata = embed_sequences(
        frame.sequence.dropna(),
        args.output,
        args.model_id,
        args.revision,
        args.batch_size,
        args.threads,
        args.local_files_only,
    )
    print(metadata)


if __name__ == "__main__":
    main()
