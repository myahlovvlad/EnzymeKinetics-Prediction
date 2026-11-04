#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Скачивание AlphaFold PDB-структур (v6) для всех уникальных UniProt ID
из датасета enzyme-kinetics (kcat / km / ki).

- Источник 1: https://alphafold.ebi.ac.uk/files/AF-{UNIPROT}-F1-model_v6.pdb
- Источник 2 (fallback для 403/404):
  https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}  ->  [0]["pdbUrl"]
- Куда: enzyme-kinetics/data/alphafold_pdb/
- Интервал: 2.0 сек между HTTP-запросами (настраивается --rate-limit)
- Подробное логирование: stdout + audit_out/download_alphafold.log
- Возобновляемо: уже валидные .pdb пропускаются; битые перекачиваются
- Манифесты: audit_out/alphafold_download_manifest.json + .csv

Запуск:
    uv run download_alphafold_pdbs.py --dry-run
    uv run download_alphafold_pdbs.py
    uv run download_alphafold_pdbs.py --rate-limit 3.0
    uv run download_alphafold_pdbs.py --only kcat km
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import signal
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR
from typing import Optional

import pandas as pd
import requests

# --------------------------------------------------------------------------- #
# Конфиг
# --------------------------------------------------------------------------- #

DATA_ROOT = DATA_DIR
OUT_DIR = DATA_ROOT / "alphafold_pdb"
LOG_DIR = OUTPUT_DIR / "audit"

TARGETS = ["kcat", "km", "ki"]
CLUSTER_SUFFIXES = ["20cluster", "40cluster", "60cluster", "80cluster", "99cluster"]

# --- AlphaFold DB ---
# v6 актуален на 2026 год; v4 снят с раздачи (403/404).
# Если EBI выпустит v7 — поменять только эту строку.
AF_VERSION = "v6"
BASE_URL = "https://alphafold.ebi.ac.uk/files"
API_URL = "https://alphafold.ebi.ac.uk/api/prediction"
FILENAME_T = "AF-{uniprot}-F1-model_{ver}.pdb"

# Контакт в User-Agent обязателен: без него файловый сервер EBI отдаёт 403.
USER_AGENT = "enzyme-kinetics-audit/1.0 (contact: myahlov@gmail.com)"

DEFAULT_RATE_LIMIT = 2.0  # сек между запросами
DEFAULT_TIMEOUT = 60.0  # сек на один HTTP-запрос
DEFAULT_RETRIES = 3  # попыток на один UniProt
MIN_PDB_BYTES = 500  # минимальный разумный размер PDB

# --------------------------------------------------------------------------- #
# Логирование
# --------------------------------------------------------------------------- #


def setup_logging(log_name: str) -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fmt = "%(asctime)s | %(levelname)-7s | %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(LOG_DIR / log_name, mode="w", encoding="utf-8"),
        ],
    )
    return logging.getLogger("af_download")


def section(logger: logging.Logger, title: str) -> None:
    logger.info("")
    logger.info("=" * 88)
    logger.info(title)
    logger.info("=" * 88)


def fmt_elapsed(seconds: float) -> str:
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m{s:02d}s"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


# --------------------------------------------------------------------------- #
# Сбор UniProt ID
# --------------------------------------------------------------------------- #


def _csv_candidates(target: str) -> list[Path]:
    base = DATA_ROOT / "Baseline" / target
    files = [
        base / f"{target}-random_trainval.csv",
        base / f"{target}-random_test.csv",
    ]
    for sfx in CLUSTER_SUFFIXES:
        files.append(base / f"{target}-random_test_sequence_{sfx}_mean_3sim.csv")
        files.append(base / f"{target}-seq_test_sequence_{sfx}.csv")
    return files


def collect_uniprot_ids(
    logger: logging.Logger, only_targets: Optional[list[str]] = None
) -> set[str]:
    ids: set[str] = set()
    targets = only_targets or TARGETS
    for target in targets:
        logger.info(f"[{target}] сбор UniProt ID")
        for f in _csv_candidates(target):
            if not f.exists():
                logger.warning(f"  пропуск (нет файла): {f.name}")
                continue
            try:
                df = pd.read_csv(f, usecols=["uniprot"], low_memory=False)
            except Exception as e:
                logger.error(f"  ! read_csv failed {f.name}: {e}")
                continue
            u = df["uniprot"].dropna().astype(str).str.strip()
            u = u[u != ""]
            n_before = len(ids)
            ids.update(u.tolist())
            logger.info(
                f"  {f.name:58s} rows={len(df):>7d} "
                f"unique_in_file={u.nunique():>7d} new={len(ids) - n_before:>7d}"
            )
    logger.info(f"Итого уникальных UniProt ID: {len(ids)}")
    return ids


# --------------------------------------------------------------------------- #
# Скачивание
# --------------------------------------------------------------------------- #


def _is_valid_pdb(text: str) -> bool:
    if len(text) < MIN_PDB_BYTES:
        return False
    return ("ATOM " in text) or ("HETATM" in text)


def _make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update(
        {
            "User-Agent": USER_AGENT,
            "Accept": "text/plain",
            "Connection": "keep-alive",
        }
    )
    return s


def _get_pdb_url_from_api(
    uniprot: str, session: requests.Session, timeout: float
) -> Optional[str]:
    """
    Возвращает актуальный pdbUrl из официального API AlphaFold DB.
    Используется как fallback для 403/404 на прямом URL.
    """
    try:
        r = session.get(f"{API_URL}/{uniprot}", timeout=timeout)
        if r.status_code == 200:
            data = r.json()
            if isinstance(data, list) and data:
                url = data[0].get("pdbUrl")
                if url:
                    return url
    except Exception:
        pass
    return None


def _atomic_save(filepath: Path, content: bytes) -> None:
    tmp = filepath.with_suffix(".pdb.part")
    tmp.write_bytes(content)
    tmp.replace(filepath)


def download_one(
    uniprot: str,
    out_dir: Path,
    session: requests.Session,
    rate_limit: float,
    timeout: float,
    retries: int,
) -> tuple[str, Optional[str]]:
    """
    Возвращает (status, detail):
        ("ok",         None)
        ("skip",       None)              уже валидный файл на месте
        ("not_found",  "404")             нет в AF DB
        ("http_error", "403")             http-ошибка
        ("invalid",    "invalid PDB (len=...)" )
        ("error",      "network: ...")
    """
    filename = FILENAME_T.format(uniprot=uniprot, ver=AF_VERSION)
    filepath = out_dir / filename

    # --- resume: уже валидный файл? ---
    if filepath.exists() and filepath.stat().st_size >= MIN_PDB_BYTES:
        try:
            head = filepath.read_text(encoding="utf-8", errors="ignore")
            if _is_valid_pdb(head):
                return "skip", None
        except Exception:
            pass
        try:
            filepath.unlink()
        except Exception:
            pass

    url = f"{BASE_URL}/{filename}"
    last_err: Optional[str] = None

    for attempt in range(1, retries + 1):
        try:
            resp = session.get(url, timeout=timeout)
        except requests.RequestException as e:
            last_err = f"network: {type(e).__name__}: {e}"
            time.sleep(rate_limit)
            continue

        if resp.status_code == 200:
            text = resp.text
            if _is_valid_pdb(text):
                _atomic_save(filepath, resp.content)
                time.sleep(rate_limit)
                return "ok", None
            last_err = f"invalid PDB (len={len(text)})"
            time.sleep(rate_limit)
            continue

        if resp.status_code in (403, 404):
            # Fallback: запросим актуальный pdbUrl из API (один раз).
            api_url = _get_pdb_url_from_api(uniprot, session, timeout)
            if api_url:
                try:
                    r2 = session.get(api_url, timeout=timeout)
                    if r2.status_code == 200 and _is_valid_pdb(r2.text):
                        _atomic_save(filepath, r2.content)
                        time.sleep(rate_limit)
                        return "ok", None
                    last_err = f"api fallback http {r2.status_code}"
                except requests.RequestException as e:
                    last_err = f"api fallback network: {type(e).__name__}"
            else:
                last_err = "api fallback: no pdbUrl"
            # 404 без валидного fallback = нет в базе, дальше не пытаемся
            if resp.status_code == 404 and last_err.startswith(
                "api fallback: no pdbUrl"
            ):
                time.sleep(rate_limit)
                return "not_found", "404"
            time.sleep(rate_limit)
            continue

        # прочие коды (429/5xx) — ретраим
        last_err = f"http {resp.status_code}"
        time.sleep(rate_limit)

    # классификация последней ошибки
    if last_err and last_err.startswith("invalid"):
        return "invalid", last_err
    if last_err and "no pdbUrl" in last_err:
        return "not_found", "404 (api)"
    return "error", last_err


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--rate-limit",
        type=float,
        default=DEFAULT_RATE_LIMIT,
        help=f"сек между запросами (по умолчанию {DEFAULT_RATE_LIMIT})",
    )
    ap.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_TIMEOUT,
        help=f"таймаут одного запроса, сек (по умолчанию {DEFAULT_TIMEOUT})",
    )
    ap.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help=f"попыток на один UniProt (по умолчанию {DEFAULT_RETRIES})",
    )
    ap.add_argument(
        "--dry-run", action="store_true", help="только посчитать, ничего не скачивать"
    )
    ap.add_argument(
        "--only",
        nargs="+",
        choices=TARGETS,
        default=None,
        help="ограничить сбор ID подмножеством target'ов",
    )
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    logger = setup_logging("download_alphafold.log")
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"AF_VERSION   = {AF_VERSION}")
    logger.info(f"USER_AGENT   = {USER_AGENT}")
    logger.info(f"OUT_DIR      = {out_dir}")
    logger.info(f"LOG          = {LOG_DIR / 'download_alphafold.log'}")
    logger.info(f"rate_limit   = {args.rate_limit} s")
    logger.info(f"timeout      = {args.timeout} s")
    logger.info(f"retries      = {args.retries}")
    logger.info(f"started      = {datetime.now().isoformat(timespec='seconds')}")

    def _sigint(signum, frame):
        raise KeyboardInterrupt()

    signal.signal(signal.SIGINT, _sigint)

    section(logger, "СБОР UNIQUE UNIPROT ID ИЗ CSV")
    try:
        ids = collect_uniprot_ids(logger, only_targets=args.only)
    except Exception:
        logger.error(traceback.format_exc())
        return 1

    if not ids:
        logger.error("Не найдено ни одного UniProt ID — проверьте DATA_ROOT и CSV.")
        return 1

    todo = sorted(ids)
    total = len(todo)

    eta_total = total * args.rate_limit
    logger.info(f"Файлов к скачиванию: {total}")
    logger.info(
        f"Ожидаемое время при rate_limit={args.rate_limit}s: ~{fmt_elapsed(eta_total)}"
    )

    if args.dry_run:
        logger.info("[DRY-RUN] скачивание не выполняется.")
        return 0

    session = _make_session()

    done_ok = 0
    done_skip = 0
    failed: list[tuple[str, str, str]] = []
    t_start = time.time()

    section(logger, "СКАЧИВАНИЕ")
    try:
        for i, u in enumerate(todo, 1):
            status, detail = download_one(
                u,
                out_dir,
                session,
                rate_limit=args.rate_limit,
                timeout=args.timeout,
                retries=args.retries,
            )
            if status == "ok":
                done_ok += 1
            elif status == "skip":
                done_skip += 1
            else:
                failed.append((u, status, detail or ""))
                logger.warning(f"  [{i}/{total}] FAIL  {u}  ({status}: {detail})")

            if i == 1 or i % 25 == 0 or i == total:
                elapsed = time.time() - t_start
                rate = i / elapsed if elapsed > 0 else 0.0
                remaining = (total - i) / rate if rate > 0 else 0.0
                logger.info(
                    f"  [{i}/{total}]  ok={done_ok}  skip={done_skip}  "
                    f"fail={len(failed)}  rate={rate:.2f}/s  "
                    f"elapsed={fmt_elapsed(elapsed)}  eta={fmt_elapsed(remaining)}"
                )
    except KeyboardInterrupt:
        logger.warning("Прервано пользователем (Ctrl-C). Сохраняю манифест...")
    finally:
        elapsed = time.time() - t_start
        section(logger, "ИТОГ")
        logger.info(f"  всего            : {total}")
        logger.info(f"  скачано новых    : {done_ok}")
        logger.info(f"  уже было (skip)  : {done_skip}")
        logger.info(f"  ошибок           : {len(failed)}")
        logger.info(f"  прошло           : {fmt_elapsed(elapsed)}")
        logger.info(f"  папка            : {out_dir}")

        manifest = {
            "started": datetime.fromtimestamp(t_start).isoformat(timespec="seconds"),
            "finished": datetime.now().isoformat(timespec="seconds"),
            "af_version": AF_VERSION,
            "user_agent": USER_AGENT,
            "out_dir": str(out_dir),
            "rate_limit_sec": args.rate_limit,
            "timeout_sec": args.timeout,
            "retries": args.retries,
            "total_unique": total,
            "ok": done_ok,
            "skip": done_skip,
            "failed_count": len(failed),
            "failed": [
                {"uniprot": u, "status": s, "detail": d} for (u, s, d) in failed
            ],
            "elapsed_sec": round(elapsed, 1),
        }
        manifest_json = LOG_DIR / "alphafold_download_manifest.json"
        with open(manifest_json, "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)
        logger.info(f"  манифест JSON    : {manifest_json}")

        manifest_csv = LOG_DIR / "alphafold_download_manifest.csv"
        with open(manifest_csv, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["uniprot", "status", "detail"])
            for u, s, d in failed:
                w.writerow([u, s, d])
        logger.info(f"  манифест CSV     : {manifest_csv}")

        ids_csv = LOG_DIR / "alphafold_uniprot_ids.csv"
        with open(ids_csv, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["uniprot"])
            for u in todo:
                w.writerow([u])
        logger.info(f"  все ID           : {ids_csv}")

    return 0 if not failed else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.stderr.write("\nInterrupted.\n")
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
