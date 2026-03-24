import os
import logging
from typing import Dict, Any
import httpx
from Bio import Entrez
from core.config import settings

logger = logging.getLogger(__name__)

class ProxyManager:
    """
    Manages global proxy settings and provides configured HTTP clients.
    """
    _instance = None
    _clients: Dict[str, httpx.Client] = {}

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ProxyManager, cls).__new__(cls)
            cls._instance.refresh_proxies()
        return cls._instance

    def refresh_proxies(self):
        """
        Updates internal state and third-party libraries based on current settings.
        """
        http_proxy = settings.http_proxy
        https_proxy = settings.https_proxy

        # Update os.environ for libraries that check it (like urllib via Biopython)
        if http_proxy:
            os.environ["HTTP_PROXY"] = http_proxy
            os.environ["http_proxy"] = http_proxy
        else:
            os.environ.pop("HTTP_PROXY", None)
            os.environ.pop("http_proxy", None)

        if https_proxy:
            os.environ["HTTPS_PROXY"] = https_proxy
            os.environ["https_proxy"] = https_proxy
        else:
            os.environ.pop("HTTPS_PROXY", None)
            os.environ.pop("https_proxy", None)

        # Update Bio.Entrez
        Entrez.email = settings.entrez_email
        # Biopython Entrez uses urllib.request.urlopen which reads from os.environ.
        # However, some versions of urllib cache the proxy handler.
        # Installing a global opener ensures all urllib calls (including Entrez) use the new proxies.
        import urllib.request
        proxies = {}
        if http_proxy:
            proxies['http'] = http_proxy
        if https_proxy:
            proxies['https'] = https_proxy
        
        opener = urllib.request.build_opener(urllib.request.ProxyHandler(proxies))
        urllib.request.install_opener(opener)
        
        # Clear existing httpx clients to force them to be recreated with new proxies
        for client in self._clients.values():
            try:
                client.close()
            except Exception:
                pass
        self._clients.clear()
        
        logger.info(f"Proxies refreshed: HTTP={http_proxy}, HTTPS={https_proxy}")

    def get_client(self, name: str = "default", timeout: float = 30.0) -> httpx.Client:
        """
        Returns an httpx.Client configured with the current proxy settings.
        """
        if name not in self._clients:
            mounts = {}
            if settings.http_proxy:
                mounts["http://"] = httpx.HTTPTransport(proxy=settings.http_proxy)
            if settings.https_proxy:
                mounts["https://"] = httpx.HTTPTransport(proxy=settings.https_proxy)

            self._clients[name] = httpx.Client(
                mounts=mounts,
                timeout=timeout
            )
        return self._clients[name]

proxy_manager = ProxyManager()
