import json
import logging
import redis
from core.config import settings

logger = logging.getLogger(__name__)

class CacheManager:
    """
    Manages connection to Redis for global caching.
    Provides a safe interface that fails gracefully if Redis is unavailable.
    """
    def __init__(self):
        self.client = None
        if settings.redis_url:
            try:
                # Using decode_responses=True to handle strings directly
                self.client = redis.from_url(settings.redis_url, decode_responses=True, socket_timeout=2.0)
                self.client.ping()
                logger.info(f"Connected to Redis at {settings.redis_url}")
            except Exception as e:
                logger.warning(f"Could not connect to Redis: {e}. Caching will be disabled.")
                self.client = None

    def get(self, key: str):
        """Retrieve a JSON-deserialized value from cache."""
        if not self.client:
            return None
        try:
            data = self.client.get(key)
            if data:
                return json.loads(data)
            return None
        except Exception as e:
            logger.error(f"Redis GET error for {key}: {e}")
            return None

    def set(self, key: str, value: any, ttl: int = 3600 * 24):
        """Store a value in cache as JSON with an optional TTL (default 24h)."""
        if not self.client:
            return
        try:
            self.client.setex(key, ttl, json.dumps(value))
        except Exception as e:
            logger.error(f"Redis SET error for {key}: {e}")

    def flush(self):
        """Clear all keys from the cache."""
        if not self.client:
            return False
        try:
            self.client.flushdb()
            logger.info("Redis cache flushed successfully.")
            return True
        except Exception as e:
            logger.error(f"Redis FLUSH error: {e}")
            return False

# Global singleton instance
cache = CacheManager()
