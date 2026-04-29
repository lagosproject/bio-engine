import os

from pydantic_settings import BaseSettings, SettingsConfigDict

from core.persistence import PersistenceManager


class Settings(BaseSettings):
    app_name: str = "Bio-Engine Sidecar"

    # Network
    host: str = "127.0.0.1"
    port: int = 8000

    # External Services
    entrez_email: str = "example@example.com"
    http_proxy: str | None = os.getenv("http_proxy") or os.getenv("HTTP_PROXY")
    https_proxy: str | None = os.getenv("https_proxy") or os.getenv("HTTPS_PROXY")

    # Paths (Defaults handled by PersistenceManager if not set)
    # We use a factory to delay PersistenceManager instantiation if possible,
    # but for simplicity we can instantiate it here as it just resolves paths.
    _persistence: PersistenceManager = PersistenceManager()

    cache_dir: str = os.getenv("BIO_CACHE_DIR", _persistence.get_cache_dir())
    jobs_dir: str = _persistence.get_jobs_dir()
    logs_dir: str = _persistence.get_logs_dir()
    uploads_dir: str = _persistence.get_uploads_dir()

    def get_cache_dir(self, assembly: str | None = None) -> str:
        if not assembly:
            return self.cache_dir
        
        # If we have an assembly, we use the persistence manager to resolve the subfolder
        # relative to the base directory of the current cache_dir
        base_cache = self.cache_dir
        norm_assembly = assembly.lower().replace("hg19", "grch37").replace("hg38", "grch38")
        return os.path.join(base_cache, norm_assembly)

    # Tracy
    tracy_path: str = os.getenv("TRACY_PATH", "tracy")

    # Samtools & Bgzip
    samtools_path: str = os.getenv("BIO_SAMTOOLS_PATH", "samtools")
    bgzip_path: str = os.getenv("BIO_BGZIP_PATH", "bgzip")

    # Logging
    log_level: str = "INFO"
    json_logs: bool = True

    # Redis Cache
    redis_url: str | None = os.getenv("BIO_REDIS_URL")

    # Database
    database_url: str = os.getenv("BIO_DATABASE_URL", f"sqlite:///{_persistence.base_dir}/bio-engine.db")

    model_config = SettingsConfigDict(env_prefix="BIO_")

settings = Settings()

