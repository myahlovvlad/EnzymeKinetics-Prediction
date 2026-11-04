"""Content-addressed artifacts; inputs and settings are part of the cache identity."""

from dataclasses import asdict, is_dataclass
import hashlib
import importlib.metadata
import json
from pathlib import Path
import tempfile

SCHEMA_VERSION = 1


def digest_file(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def content_key(value):
    if is_dataclass(value):
        value = asdict(value)
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def package_version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


def artifact(
    backend, inputs, settings, values, units, *, limitations=(), status="complete"
):
    if set(values) != set(units):
        raise ValueError("Every physical result requires an explicit unit")
    return {
        "schema_version": SCHEMA_VERSION,
        "backend": backend,
        "status": status,
        "inputs": inputs,
        "settings": settings,
        "cache_key": content_key(
            {
                "schema": SCHEMA_VERSION,
                "backend": backend,
                "inputs": inputs,
                "settings": settings,
            }
        ),
        "values": values,
        "units": units,
        "limitations": list(limitations),
    }
