"""
OpenCRAVAT management service.

Wraps the `oc` CLI for module listing, install/uninstall, disk usage reporting,
and ad-hoc variant annotation. Install tasks run in background threads so the
API returns immediately; poll /opencravat/tasks/{task_id} for progress.
"""

import logging
import re
import shutil
import subprocess
import threading
import uuid
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# In-memory install/uninstall task registry. Survives for the process lifetime.
_tasks: dict[str, dict[str, Any]] = {}

_SIZE_PATTERN = re.compile(r"([\d.]+)\s*(B|KB|MB|GB|TB)", re.IGNORECASE)
_SIZE_UNITS = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _run(
    *args: str,
    oc_path: str = "oc",
    timeout: int = 60,
    stdin_input: str | None = None,
) -> tuple[int, str, str]:
    try:
        result = subprocess.run(
            [oc_path, *args],
            capture_output=True,
            text=True,
            timeout=timeout,
            input=stdin_input,
        )
        return result.returncode, result.stdout.strip(), result.stderr.strip()
    except FileNotFoundError:
        return (
            -1,
            "",
            f"OpenCRAVAT binary not found at '{oc_path}'. "
            "Install with: pip install open-cravat && oc install-base",
        )
    except subprocess.TimeoutExpired:
        return -1, "", f"Command timed out after {timeout}s"


def _parse_size_bytes(size_str: str) -> int | None:
    m = _SIZE_PATTERN.search(size_str)
    if not m:
        return None
    return int(float(m.group(1)) * _SIZE_UNITS.get(m.group(2).upper(), 1))


def _dir_size_bytes(path: str) -> int:
    total = 0
    for f in Path(path).rglob("*"):
        if f.is_file():
            try:
                total += f.stat().st_size
            except OSError:
                pass
    return total


def _get_data_dir(oc_path: str) -> str | None:
    rc, out, _ = _run("config", "md", oc_path=oc_path)
    if rc == 0 and out:
        return out
    default = Path.home() / ".open-cravat"
    return str(default) if default.exists() else None


def _parse_module_table(stdout: str, default_installed: bool | None = None) -> list[dict[str, Any]]:
    """
    Parse the tabular output of `oc module ls` or `oc module ls -a`.

    OC output columns (order may vary by version):
        name  version  type  data_size  [yes/no]  title
    """
    known_types = {"annotator", "converter", "postaggregator", "reporter", "mapper"}
    modules: list[dict[str, Any]] = []

    for line in stdout.splitlines():
        line = line.strip()
        if not line or re.match(r"^[-=\s]*$", line):
            continue
        # Skip header lines that contain column labels
        if re.match(r"(?i)^\s*(name|module)\s", line):
            continue

        parts = line.split()
        if not parts:
            continue

        name = parts[0]
        mod: dict[str, Any] = {
            "name": name,
            "version": None,
            "type": None,
            "size_bytes": None,
            "installed": default_installed,
            "title": None,
        }

        # Version: first token that looks like a version number
        for p in parts[1:]:
            if re.match(r"^\d[\d.]{1,}", p):
                mod["version"] = p
                break

        # Type: known module type keyword
        for p in parts[1:]:
            if p.lower() in known_types:
                mod["type"] = p.lower()
                break

        # Installed: explicit yes/no column
        for p in parts[1:]:
            if p.lower() == "yes":
                mod["installed"] = True
                break
            if p.lower() == "no":
                mod["installed"] = False
                break

        # Size: number immediately followed by a unit token
        for i, p in enumerate(parts[:-1]):
            if re.match(r"^[\d.]+$", p) and parts[i + 1].upper() in _SIZE_UNITS:
                mod["size_bytes"] = _parse_size_bytes(f"{p} {parts[i + 1]}")
                break

        modules.append(mod)

    return modules


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_status(oc_path: str = "oc") -> dict[str, Any]:
    rc, version_out, err = _run("--version", oc_path=oc_path)
    if rc != 0:
        return {"installed": False, "error": err}

    data_dir = _get_data_dir(oc_path)
    disk_used_bytes: int | None = None
    disk_free_bytes: int | None = None

    if data_dir and Path(data_dir).exists():
        disk_used_bytes = _dir_size_bytes(data_dir)
        disk_free_bytes = shutil.disk_usage(data_dir).free

    return {
        "installed": True,
        "version": version_out,
        "data_dir": data_dir,
        "disk_used_bytes": disk_used_bytes,
        "disk_free_bytes": disk_free_bytes,
        "error": None,
    }


def list_installed(oc_path: str = "oc") -> list[dict[str, Any]]:
    rc, stdout, stderr = _run("module", "ls", oc_path=oc_path, timeout=30)
    if rc != 0:
        logger.error(f"oc module ls failed: {stderr}")
        return []
    return _parse_module_table(stdout, default_installed=True)


def list_store(oc_path: str = "oc") -> list[dict[str, Any]]:
    """List all modules available in the OC store (installed + available)."""
    rc, stdout, stderr = _run("module", "ls", "-a", oc_path=oc_path, timeout=60)
    if rc != 0:
        logger.error(f"oc module ls -a failed: {stderr}")
        return []
    return _parse_module_table(stdout)


def get_task(task_id: str) -> dict[str, Any] | None:
    return _tasks.get(task_id)


def _do_install(task_id: str, module_name: str, oc_path: str) -> None:
    _tasks[task_id]["status"] = "running"
    rc, _, stderr = _run(
        "module", "install", "-y", module_name,
        oc_path=oc_path,
        timeout=3600,   # large installs (gnomAD ~60 GB) can take a long time
        stdin_input="y\n",
    )
    if rc == 0:
        _tasks[task_id]["status"] = "completed"
        logger.info(f"OpenCRAVAT module '{module_name}' installed successfully")
    else:
        _tasks[task_id]["status"] = "failed"
        _tasks[task_id]["error"] = stderr
        logger.error(f"OpenCRAVAT install '{module_name}' failed: {stderr[:300]}")


def install_module(module_name: str, oc_path: str = "oc") -> str:
    """Start a background install. Returns task_id to poll for status."""
    task_id = uuid.uuid4().hex[:8]
    _tasks[task_id] = {
        "task_id": task_id,
        "module": module_name,
        "status": "pending",
        "error": None,
    }
    threading.Thread(
        target=_do_install, args=(task_id, module_name, oc_path), daemon=True
    ).start()
    return task_id


def uninstall_module(module_name: str, oc_path: str = "oc") -> dict[str, Any]:
    rc, _, stderr = _run("module", "uninstall", "-y", module_name, oc_path=oc_path, timeout=120)
    if rc == 0:
        logger.info(f"OpenCRAVAT module '{module_name}' uninstalled")
        return {"success": True, "error": None}
    return {"success": False, "error": stderr}
