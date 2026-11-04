"""Append-only, process-safe experiment event journal with immutable registrations."""

from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import time
import uuid

from enzymekinetics_prediction.physics.artifacts import content_key, write_json


class Journal:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)

    def event(self, experiment, event, **details):
        record = {
            "utc": datetime.now(timezone.utc).isoformat(),
            "experiment": experiment,
            "event": event,
            **details,
        }
        serialized = json.dumps(record, ensure_ascii=False, allow_nan=False)
        with (self.directory / "events.jsonl").open("a", encoding="utf-8") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            handle.write(serialized + "\n")
            handle.flush()
            import os

            os.fsync(handle.fileno())
            fcntl.flock(handle, fcntl.LOCK_UN)
        return record

    def register(self, hypotheses):
        path = self.directory / "registration.json"
        if path.exists():
            if json.loads(path.read_text()) != hypotheses:
                raise ValueError(
                    "A registration cannot be changed after experiments start"
                )
            return
        # Exclusive creation prevents independent runners from changing registration.
        with path.open("x") as handle:
            json.dump(hypotheses, handle, ensure_ascii=False, indent=2, allow_nan=False)
        self.event("registration", "registered", sha256=content_key(hypotheses))

    @contextmanager
    def run(self, experiment, parameters, inputs=None):
        run_id = str(uuid.uuid4())
        start = time.monotonic()
        self.event(
            experiment,
            "started",
            run_id=run_id,
            parameters=parameters,
            inputs=inputs or {},
        )
        try:
            yield run_id
        except BaseException as exc:
            self.event(
                experiment,
                "failed",
                run_id=run_id,
                seconds=time.monotonic() - start,
                error_type=type(exc).__name__,
                error=str(exc),
            )
            raise
        else:
            self.event(
                experiment, "completed", run_id=run_id, seconds=time.monotonic() - start
            )

    def result(self, experiment, result):
        history = self.directory / "results"
        history.mkdir(exist_ok=True)
        latest = self.directory / (experiment + ".json")
        if latest.exists():
            previous = json.loads(latest.read_text())
            write_json(
                history / (experiment + "-" + content_key(previous) + ".json"), previous
            )
        immutable = history / (experiment + "-" + content_key(result) + ".json")
        write_json(immutable, result)
        write_json(self.directory / (experiment + ".json"), result)
        self.event(
            experiment,
            "result",
            result_sha256=content_key(result),
            artifact=str(immutable.relative_to(self.directory)),
        )
