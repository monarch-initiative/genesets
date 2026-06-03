from __future__ import annotations

import datetime as dt
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def repo_root(start: Path | None = None) -> Path:
    path = (start or Path.cwd()).resolve()
    for candidate in (path, *path.parents):
        if (candidate / "Cargo.toml").exists() and (candidate / "src").exists():
            return candidate
    return path


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)


def resolve_path(base: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def run_command(command: list[str], cwd: Path | None = None) -> dict[str, Any]:
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    elapsed = round(time.perf_counter() - started, 3)
    if completed.returncode != 0:
        print(completed.stdout, file=sys.stdout)
        print(completed.stderr, file=sys.stderr)
        raise SystemExit(
            f"command failed with exit code {completed.returncode}: {' '.join(command)}"
        )
    return {
        "command": command,
        "runtime_seconds": elapsed,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def duckdb_json(sql: str) -> list[dict[str, Any]] | None:
    if not shutil.which("duckdb"):
        return None
    completed = subprocess.run(
        ["duckdb", "-json", "-c", sql],
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(completed.stdout or "[]")


def duckdb_python_query(path: Path, sql: str) -> list[tuple[Any, ...]] | None:
    try:
        import duckdb  # type: ignore
    except ImportError:
        return None

    con = duckdb.connect()
    return con.execute(sql, [str(path)]).fetchall()


def which(command: str) -> str | None:
    return shutil.which(command)


def version_output(command: list[str]) -> str | None:
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    output = completed.stdout.strip() or completed.stderr.strip()
    return output or None
