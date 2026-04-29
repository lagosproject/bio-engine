import logging
"""
Persistence Management Module
=============================

This module provides the `PersistenceManager` class, responsible for 
resolving and creating the necessary local directories (logs, jobs, cache) 
based on the host operating system conventions.
"""

import os
import pathlib
import sys

logger = logging.getLogger(__name__)

class PersistenceManager:
    """
    Manages local storage paths for the application, ensuring that 
    required directories exist across different operating systems.
    """
    def __init__(self, app_name: str = "ps-analyzer"):
        self.app_name = app_name
        self.base_dir = self._resolve_base_dir()
        self._ensure_dirs()

    def _resolve_base_dir(self) -> pathlib.Path:
        # 1. Environment variable override (highest priority)
        if env_dir := os.environ.get("BIO_DATA_DIR"):
            return pathlib.Path(env_dir)

        # 2. Check for portable mode
        # We consider it portable if BIO_PORTABLE is set OR if a .portable file exists next to the exe
        exe_dir = self._get_exe_dir()
        is_portable = os.environ.get("BIO_PORTABLE") == "1"
        
        if not is_portable:
            # Check for marker file in the executable directory or its parent (if in 'binaries' subfolder)
            if (exe_dir / ".portable").exists() or (exe_dir / "portable").exists():
                is_portable = True
            elif exe_dir.name == "binaries" and ((exe_dir.parent / ".portable").exists() or (exe_dir.parent / "portable").exists()):
                is_portable = True

        if is_portable:
            if exe_dir.name == "binaries":
                return exe_dir.parent / "data"
            return exe_dir / "data"

        # 3. AppImage specific handling
        if sys.platform == "linux" and "APPIMAGE" in os.environ:
            appimage_dir = pathlib.Path(os.environ["APPIMAGE"]).parent
            # Only use it if it's writable, otherwise fall back to standard
            if os.access(appimage_dir, os.W_OK):
                return appimage_dir / f"{self.app_name}-data"

        # 4. Standard OS-specific locations
        home = pathlib.Path.home()
        if sys.platform == "win32":
            return home / "AppData" / "Roaming" / self.app_name
        elif sys.platform == "darwin":
            return home / "Library" / "Application Support" / self.app_name
        else:  # Linux/Unix
            return home / ".local" / "share" / self.app_name

    def _get_exe_dir(self) -> pathlib.Path:
        if getattr(sys, 'frozen', False):
            return pathlib.Path(sys.executable).parent
        return pathlib.Path(__file__).parent.parent

    def _ensure_dirs(self):
        dirs = [
            self.get_jobs_dir(),
            self.get_logs_dir(),
            self.get_cache_dir(),
            self.get_uploads_dir()
        ]
        for d in dirs:
            try:
                os.makedirs(d, exist_ok=True)
            except OSError as e:
                logger.error(f"Could not create directory {d}: {e}")
                # Fallback to local 'data' directory relative to the current file
                fallback = pathlib.Path(__file__).parent.parent / "data" / os.path.basename(d)
                try:
                    os.makedirs(fallback, exist_ok=True)
                    logger.warning(f"Falling back to {fallback}")
                except Exception as ex:
                    logger.critical(f"Failed to create even fallback directory: {ex}")

    def get_jobs_dir(self) -> str:
        return str(self.base_dir / "jobs")

    def get_logs_dir(self) -> str:
        return str(self.base_dir / "logs")

    def get_cache_dir(self, assembly: str | None = None) -> str:
        base_cache = str(self.base_dir / "ncbi_cache")
        if assembly:
            # Normalize assembly name (e.g. hg38, GRCh38 -> grch38)
            norm_assembly = assembly.lower().replace("hg19", "grch37").replace("hg38", "grch38")
            return os.path.join(base_cache, norm_assembly)
        return base_cache

    def get_uploads_dir(self) -> str:
        return str(self.base_dir / "uploads")

    def get_log_file(self, filename: str = "bio-engine.log") -> str:
        return os.path.join(self.get_logs_dir(), filename)
