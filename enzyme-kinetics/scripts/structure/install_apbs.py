"""Install the pinned official APBS Linux release locally, verifying its SHA-256."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import tempfile
import urllib.request
import zipfile

from enzymekinetics_prediction.paths import VENDOR_DIR
from enzymekinetics_prediction.physics.artifacts import digest_file, write_json

VERSION = "3.5.0"
URL = "https://github.com/Electrostatics/apbs/releases/download/v3.5.0/APBS-3.5.0.Linux.zip"
SHA256 = "b92f6a1bc05c9d4278d85a68072ab6b1ce5aed8d4bb53870746104635b923ca3"
EXECUTABLE = Path("APBS-3.5.0.Linux/bin/apbs")


def install(destination):
    destination = Path(destination)
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "AMD64"}:
        raise ValueError(
            "Pinned package is Linux x86_64 only; use APBS_EXECUTABLE on other platforms"
        )
    manifest = destination / "installation.json"
    if destination.exists() and any(destination.iterdir()):
        if manifest.is_file():
            old = json.loads(manifest.read_text())
            if old.get("archive_sha256") == SHA256 and digest_file(
                destination / EXECUTABLE
            ) == old.get("executable_sha256"):
                return old
        raise ValueError("Refusing to overwrite an unverified existing installation")
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="enzyme-apbs-") as temporary:
        archive = Path(temporary) / "apbs.zip"
        digest, size = hashlib.sha256(), 0
        with (
            urllib.request.urlopen(URL, timeout=120) as response,
            archive.open("wb") as target,
        ):
            while chunk := response.read(1024 * 1024):
                size += len(chunk)
                if size > 80_000_000:
                    raise ValueError("Download exceeds expected release size bound")
                digest.update(chunk)
                target.write(chunk)
        if digest.hexdigest() != SHA256:
            raise ValueError("Official APBS archive checksum mismatch")
        with zipfile.ZipFile(archive) as bundle:
            if sum(item.file_size for item in bundle.infolist()) > 1_000_000_000:
                raise ValueError("Unpacked archive exceeds limit")
            for item in bundle.infolist():
                path = (destination / item.filename).resolve()
                if not path.is_relative_to(destination.resolve()):
                    raise ValueError("Unsafe archive path")
            bundle.extractall(destination)
            for item in bundle.infolist():
                mode = item.external_attr >> 16
                if mode:
                    (destination / item.filename).chmod(mode & 0o777)
    result = {
        "version": VERSION,
        "source": URL,
        "archive_sha256": SHA256,
        "executable": str(EXECUTABLE),
        "executable_sha256": digest_file(destination / EXECUTABLE),
        "license": "Retained in official extracted distribution; see share/apbs",
    }
    write_json(manifest, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=VENDOR_DIR / "apbs")
    args = parser.parse_args()
    print(json.dumps(install(args.destination)))


if __name__ == "__main__":
    main()
