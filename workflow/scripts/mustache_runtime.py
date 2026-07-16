"""Socket-free execution support for the pinned Mustache 1.3.3 caller.

Mustache imports ``multiprocessing.Manager`` and ``multiprocessing.Process``
directly into :mod:`mustache.mustache`.  A Manager starts a local socket server,
which is unnecessary for Mustache's block workers and is unavailable on some
secured compute nodes.  The classes below provide the tiny API Mustache uses,
while keeping its existing batching and worker-count logic intact.
"""

import math
import sys
import threading
import traceback
from contextvars import ContextVar, Token
from pathlib import Path
from types import ModuleType
from typing import Any, Callable, Iterable, Iterator, Sequence


MUSTACHE_COLUMNS = (
    "BIN1_CHR",
    "BIN1_START",
    "BIN1_END",
    "BIN2_CHROMOSOME",
    "BIN2_START",
    "BIN2_END",
    "FDR",
    "DETECTION_SCALE",
)
MUSTACHE_HEADER = "\t".join(MUSTACHE_COLUMNS) + "\n"
EXECUTION_BACKEND = "manager-free-threads"


class MustacheWorkerError(RuntimeError):
    """A Mustache thread failed while processing one contact-map block."""


class _SharedList:
    """The synchronized list subset used by Mustache's block workers."""

    def __init__(self, values: Iterable[Any] = ()) -> None:
        self._values = list(values)
        self._lock = threading.Lock()

    def append(self, value: Any) -> None:
        with self._lock:
            self._values.append(value)

    def __len__(self) -> int:
        with self._lock:
            return len(self._values)

    def __getitem__(self, index: int | slice) -> Any:
        with self._lock:
            return self._values[index]

    def __iter__(self) -> Iterator[Any]:
        with self._lock:
            snapshot = tuple(self._values)
        return iter(snapshot)


class _ThreadGroup:
    """Track an invocation's workers so exceptions cannot orphan live threads."""

    def __init__(self) -> None:
        self.processes: list["ThreadProcess"] = []

    def add(self, process: "ThreadProcess") -> None:
        self.processes.append(process)

    def drain(self) -> list[MustacheWorkerError]:
        errors = []
        for process in self.processes:
            if not process.started:
                continue
            try:
                process.join()
            except MustacheWorkerError as error:
                errors.append(error)
        return errors


_ACTIVE_GROUP: ContextVar[_ThreadGroup | None] = ContextVar(
    "mustache_thread_group", default=None
)


class ThreadManager:
    """Context manager compatible with the API Mustache uses from Manager."""

    def __init__(self) -> None:
        self._group = _ThreadGroup()
        self._token: Token[_ThreadGroup | None] | None = None

    def __enter__(self) -> "ThreadManager":
        self._token = _ACTIVE_GROUP.set(self._group)
        return self

    def __exit__(self, exc_type, exc, exc_tb) -> bool:  # type: ignore[no-untyped-def]
        try:
            errors = self._group.drain()
        finally:
            if self._token is not None:
                _ACTIVE_GROUP.reset(self._token)
                self._token = None
        if exc_type is None and errors:
            raise errors[0]
        return False

    def list(self, values: Iterable[Any] = ()) -> _SharedList:
        return _SharedList(values)


class ThreadProcess:
    """Small ``multiprocessing.Process`` analogue backed by one Python thread."""

    def __init__(
        self,
        group: Any = None,
        target: Callable[..., Any] | None = None,
        name: str | None = None,
        args: Sequence[Any] = (),
        kwargs: dict[str, Any] | None = None,
        *,
        daemon: bool | None = None,
    ) -> None:
        if group is not None:
            raise ValueError("group argument must be None")
        if target is None:
            raise ValueError("a worker target is required")
        self._target = target
        self._args = tuple(args)
        self._kwargs = dict(kwargs or {})
        self._failure: tuple[BaseException, str] | None = None
        self._started = False
        self._thread = threading.Thread(
            target=self._run,
            name=name,
            daemon=bool(daemon) if daemon is not None else False,
        )
        active_group = _ACTIVE_GROUP.get()
        if active_group is not None:
            active_group.add(self)

    def _run(self) -> None:
        try:
            self._target(*self._args, **self._kwargs)
        except BaseException as error:  # propagate failures at join, like a process exit
            self._failure = (error, traceback.format_exc())

    @property
    def name(self) -> str:
        return self._thread.name

    @property
    def started(self) -> bool:
        return self._started

    @property
    def exitcode(self) -> int | None:
        if self._thread.ident is None or self._thread.is_alive():
            return None
        return 1 if self._failure is not None else 0

    def start(self) -> None:
        self._thread.start()
        self._started = True

    def is_alive(self) -> bool:
        return self._thread.is_alive()

    def join(self, timeout: float | None = None) -> None:
        self._thread.join(timeout)
        if self._thread.is_alive() or self._failure is None:
            return
        error, worker_traceback = self._failure
        raise MustacheWorkerError(
            f"Mustache worker {self.name!r} failed:\n{worker_traceback.rstrip()}"
        ) from error


def run_mustache_threaded(
    cli_args: Sequence[str], module: ModuleType | Any | None = None
) -> None:
    """Run Mustache's normal CLI after patching only its concurrency primitives."""
    if module is None:
        from mustache import mustache as module  # type: ignore[no-redef]

    if not hasattr(module, "Manager") or not hasattr(module, "Process"):
        raise RuntimeError(
            "Installed Mustache is incompatible: Manager/Process imports are missing"
        )

    previous_manager = module.Manager
    previous_process = module.Process
    previous_argv = sys.argv
    module.Manager = ThreadManager
    module.Process = ThreadProcess
    sys.argv = ["mustache", *map(str, cli_args)]
    try:
        return_code = module.main()
    finally:
        sys.argv = previous_argv
        module.Manager = previous_manager
        module.Process = previous_process

    if return_code not in (None, 0):
        raise RuntimeError(f"Mustache exited with status {return_code!r}")


def _chromosome_key(chromosome: str) -> tuple[int, int | str, str]:
    label = chromosome[3:] if chromosome.lower().startswith("chr") else chromosome
    upper = label.upper()
    if label.isdigit():
        return (0, int(label), chromosome)
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if upper in special:
        return (0, special[upper], chromosome)
    return (1, upper, chromosome)


def validate_and_sort_output(path: str | Path, resolution: int) -> int:
    """Validate Mustache's BEDPE-like table and rewrite it in genomic order."""
    if resolution <= 0:
        raise ValueError("Mustache resolution must be positive")

    output = Path(path)
    if not output.is_file():
        raise RuntimeError(f"Mustache did not create its output: {output}")
    lines = output.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise RuntimeError("Mustache created an empty output file")
    if tuple(lines[0].split("\t")) != MUSTACHE_COLUMNS:
        raise RuntimeError("Mustache output has an unexpected header")

    records: list[tuple[tuple[Any, ...], tuple[str, ...]]] = []
    for line_number, line in enumerate(lines[1:], start=2):
        if not line:
            raise RuntimeError(f"Mustache output has a blank row at line {line_number}")
        fields = tuple(line.split("\t"))
        if len(fields) != len(MUSTACHE_COLUMNS):
            raise RuntimeError(
                f"Mustache output line {line_number} has {len(fields)} columns; "
                f"expected {len(MUSTACHE_COLUMNS)}"
            )
        chrom1, chrom2 = fields[0], fields[3]
        if not chrom1 or not chrom2 or chrom1 != chrom2:
            raise RuntimeError(
                f"Mustache output line {line_number} is not a valid cis loop"
            )
        try:
            start1, end1 = int(fields[1]), int(fields[2])
            start2, end2 = int(fields[4]), int(fields[5])
            fdr, scale = float(fields[6]), float(fields[7])
        except ValueError as error:
            raise RuntimeError(
                f"Mustache output line {line_number} contains a non-numeric value"
            ) from error
        if (
            start1 < 0
            or start2 < 0
            or end1 - start1 != resolution
            or end2 - start2 != resolution
            or start1 % resolution
            or start2 % resolution
        ):
            raise RuntimeError(
                f"Mustache output line {line_number} is off the {resolution}-bp grid"
            )
        if not math.isfinite(fdr) or not 0.0 <= fdr <= 1.0:
            raise RuntimeError(
                f"Mustache output line {line_number} has an invalid FDR"
            )
        if not math.isfinite(scale) or scale <= 0.0:
            raise RuntimeError(
                f"Mustache output line {line_number} has an invalid detection scale"
            )
        key = (
            _chromosome_key(chrom1),
            start1,
            end1,
            _chromosome_key(chrom2),
            start2,
            end2,
            fdr,
            scale,
            fields,
        )
        records.append((key, fields))

    records.sort(key=lambda item: item[0])
    canonical = [MUSTACHE_HEADER]
    canonical.extend("\t".join(fields) + "\n" for _, fields in records)
    output.write_text("".join(canonical), encoding="utf-8")
    return len(records)
