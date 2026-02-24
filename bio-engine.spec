# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all
import sys
import os

# Get the environment prefix dynamically
conda_prefix = sys.prefix

datas = []
binaries = []

import shutil

# Tracy will be handled by Tauri, not bundled here


# Optionally include these shared libraries only if they exist (Linux/Conda specific fixes)
for lib in ['libcrypto.so.3', 'libssl.so.3']:
    lib_path = os.path.join(conda_prefix, 'lib', lib)
    if os.path.exists(lib_path):
        binaries.append((lib_path, '.'))

hiddenimports = [
    'fastapi', 
    'uvicorn', 
    'uvicorn.logging',
    'uvicorn.loops',
    'uvicorn.loops.auto',
    'uvicorn.protocols',
    'uvicorn.protocols.http',
    'uvicorn.protocols.http.auto',
    'uvicorn.protocols.http.httptools_impl', 
    'uvicorn.protocols.http.h11_impl', 
    'uvicorn.protocols.websockets',
    'uvicorn.protocols.websockets.auto',
    'uvicorn.protocols.websockets.websockets_impl', 
    'uvicorn.lifespan',
    'uvicorn.lifespan.on',
    'uvicorn.lifespan.off',
    'uvicorn.config',
    'uvicorn.main',
    'uvicorn.server',
    'uvicorn._types',
    'h11',
    'h11._connection',
    'h11._events',
    'h11._state',
    'h11._readers',
    'h11._writers',
    'h11._util',
    'h11._abnf',
    'h11._headers',
    'h11._receivebuffer',
    'click',
    'anyio',
    'anyio._backends',
    'anyio._backends._asyncio',
    'sniffio',
    'starlette',
    'starlette.routing',
    'starlette.responses',
    'starlette.requests',
    'starlette.middleware',
    'starlette.middleware.cors',
    'psycopg2',
    'pkg_resources'
]

# Collect all dependencies
packages_to_collect = [
    'fastapi', 'uvicorn', 'h11', 'click', 'anyio', 'sniffio', 'starlette',
    'Bio', 'hgvs', 'ometa', 'parsley', 'terml', 
    'psycopg2', 'bioutils', 'pkg_resources'
]

for package in packages_to_collect:
    try:
        tmp_ret = collect_all(package)
        datas += tmp_ret[0]
        binaries += tmp_ret[1]
        hiddenimports += tmp_ret[2]
    except Exception as e:
        print(f"WARNING: Could not collect {package}: {e}")

a = Analysis(
    ['main.py'],
    pathex=['.', os.path.join(conda_prefix, 'lib')],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='bio-engine',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
