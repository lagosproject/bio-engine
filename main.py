"""
Bio-Engine Sidecar Entry Point
==============================

This module is the entry point for the PS Analyzer Bio-Engine sidecar.
It initializes the FastAPI application, sets up logging, CORS, and
global exception handlers, and configures Uvicorn for serving the API.
"""

import argparse
import logging
import multiprocessing
import os
import signal
import sys

import uvicorn
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from api.routes import router
from core.config import settings
from core.exceptions import BioEngineError
from core.logging import setup_logging

# Handle PyInstaller one-file mode: Add temp dir to PATH so tracy is found
if getattr(sys, 'frozen', False):
    bundle_dir = sys._MEIPASS
    os.environ["PATH"] = bundle_dir + os.pathsep + os.environ["PATH"]

# Initialize Logging
setup_logging()
logger = logging.getLogger(__name__)

if getattr(sys, 'frozen', False):
    bundle_dir = sys._MEIPASS
    exe_dir = os.path.dirname(sys.executable)
    logger.info(f"Running in frozen mode.")
    logger.info(f"Bundle dir (_MEIPASS): {bundle_dir}")
    logger.info(f"Executable dir: {exe_dir}")
    
    # Add bundle_dir to PATH (for internal PyInstaller files)
    os.environ["PATH"] = bundle_dir + os.pathsep + os.environ["PATH"]
    
    # Try to find tracy and bgzip in the executable's directory
    # They might have the ps-analyzer- prefix and the target triple
    target_triple = "x86_64-pc-windows-msvc" if sys.platform == "win32" else "x86_64-unknown-linux-gnu"
    
    for tool_name, setting_attr in [("tracy", "tracy_path"), ("bgzip", "bgzip_path"), ("samtools", "samtools_path")]:
        potential_names = [
            f"ps-analyzer-{tool_name}-{target_triple}.exe" if sys.platform == "win32" else f"ps-analyzer-{tool_name}-{target_triple}",
            f"ps-analyzer-{tool_name}.exe" if sys.platform == "win32" else f"ps-analyzer-{tool_name}",
            f"{tool_name}.exe" if sys.platform == "win32" else tool_name
        ]
        
        for name in potential_names:
            path = os.path.join(exe_dir, name)
            if os.path.exists(path):
                logger.info(f"Auto-detected {tool_name} sidecar at: {path}")
                setattr(settings, setting_attr, path)
                break

app = FastAPI(title="Bio-Engine Sidecar")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.exception_handler(BioEngineError)
async def bio_engine_exception_handler(request: Request, exc: BioEngineError):
    """
    Global exception handler for BioEngineError.
    
    Transforms internal application exceptions into a standardized JSON response 
    with a 500 status code, exposing the exception type, message, and context.
    """
    return JSONResponse(
        status_code=500,
        content={
            "type": type(exc).__name__,
            "message": exc.message,
            "context": exc.context
        }
    )

# Include the API router
app.include_router(router)

def signal_handler(sig, frame):
    """
    Handles OS signals for graceful shutdown.
    
    Args:
        sig: The signal number.
        frame: The current stack frame.
    """
    logger.info(f"Received signal {sig}, shutting down...")
    sys.exit(0)

if __name__ == "__main__":
    multiprocessing.freeze_support()

    parser = argparse.ArgumentParser(description="Bio-Engine Sidecar")
    parser.add_argument("--tracy-path", type=str, help="Path to the tracy binary")
    parser.add_argument("--samtools-path", type=str, help="Path to the samtools binary")
    parser.add_argument("--bgzip-path", type=str, help="Path to the bgzip binary")
    parser.add_argument("--resource-path", type=str, help="Path to the application resources (where DLLs are located)")
    args, unknown = parser.parse_known_args()

    def clean_win_path(p):
        if p and p.startswith('\\\\?\\'):
            return p[4:]
        return p

    if args.resource_path:
        resource_path = clean_win_path(args.resource_path)
        logger.info(f"Adding resource path to PATH: {resource_path}")
        
        # Add primary resource path
        if resource_path not in os.environ["PATH"]:
            os.environ["PATH"] = resource_path + os.pathsep + os.environ["PATH"]
            
        # Also add "binaries" subfolder if it exists (common in Tauri resource structure)
        binaries_path = os.path.join(resource_path, "binaries")
        if os.path.exists(binaries_path) and binaries_path not in os.environ["PATH"]:
            logger.info(f"Adding binaries subdirectory to PATH: {binaries_path}")
            os.environ["PATH"] = binaries_path + os.pathsep + os.environ["PATH"]

    # Also add the executable's directory to PATH for sidecar DLLs
    if getattr(sys, 'frozen', False):
        exe_dir = os.path.dirname(sys.executable)
        if exe_dir not in os.environ["PATH"]:
            os.environ["PATH"] = exe_dir + os.pathsep + os.environ["PATH"]
            logger.info(f"Added executable dir to PATH: {exe_dir}")

    # Log the final path on Windows for debugging if needed
    if sys.platform == "win32":
        logger.debug(f"Effective system PATH: {os.environ['PATH']}")

    logger.info(f"Received CLI arguments: {sys.argv}")
    logger.info(f"Unknown arguments: {unknown}")
    logger.info(f"Environment BIO_SAMTOOLS_PATH: {os.getenv('BIO_SAMTOOLS_PATH')}")
    logger.info(f"Environment BIO_BGZIP_PATH: {os.getenv('BIO_BGZIP_PATH')}")
    logger.info(f"Environment TRACY_PATH: {os.getenv('TRACY_PATH')}")

    if args.tracy_path:
        logger.info(f"Overriding tracy_path from CLI: {args.tracy_path}")
        settings.tracy_path = args.tracy_path

    if args.samtools_path:
        logger.info(f"Overriding samtools_path from CLI: {args.samtools_path}")
        settings.samtools_path = args.samtools_path

    if args.bgzip_path:
        logger.info(f"Overriding bgzip_path from CLI: {args.bgzip_path}")
        settings.bgzip_path = args.bgzip_path

    # Register signal handlers for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # By default, uvicorn uses its own logging config.
    # To include uvicorn logs in our file, we pass log_config=None to use our root logger setup.
    try:
        uvicorn.run(app, host=settings.host, port=settings.port, log_config=None)
    except Exception as e:
        logger.error(f"Failed to start server: {e}")
        sys.exit(1)
