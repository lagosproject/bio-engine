# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2026-06-09
### Added
- Add local coordinate mapper for RefSeqGene assemblies
- Genericizing local annotator framework: fix OpenCRAVAT assembly casing, add fallback coordinate parsing from SPDI/HGVSG, increase Ensembl client timeout to 120s, and implement structured JSON reporter output parsing
- Integrate OpenCRAVAT module management and version endpoint
- Opencravat

### Changed
- Optimize coordinate recoder: add local parsing for NC_ genomic HGVS variants to bypass Ensembl API
- Optimize Ensembl variant recoding performance: reduce client timeout to 15s and increase chunk size to 100
- Programmatically opt out of OpenCRAVAT metrics/funding email prompts on startup to prevent hangs in docker containers
- Make VEP/OpenCRAVAT annotation cache keys mode-aware to prevent stale results when switching engines

### Fixed
- Fix OpenCRAVAT annotation parsing, integrate ClinVar, optimize HGVS coordinates, and add unit tests

## [v0.7.9] - 2026-05-28
## [v0.7.8] - 2026-05-27
### Changed
- Create release before uploading assets in pipeline

## [v0.7.7] - 2026-05-26
## [v0.7.6] - 2026-05-26
### Added
- Better cache system
- Retrieval by assembly
- Bad trimming errors
- Hotspost variants v1
- Redis cache on web
- Revied flags
- Download tracy architecture
- Update dockerfile
- Browser dialogs
- Add support for Docker server/intranet deployment

### Changed
- Refactor
- Docker images changed

### Fixed
- Del ins annotation
- Missing import

## [v0.7.5] - 2026-04-22
### Added
- Improve portable mode and path resolution

## [v0.7.4] - 2026-04-13
### Fixed
- Order on result

## [v0.7.3] - 2026-04-07
### Fixed
- Fix release assets and standardize naming for Tauri

## [v0.7.2] - 2026-04-06
### Added
- Implement aggressive DLL discovery for Windows
- Add DLL diagnostic and bundling for Windows test
- Add isolated bgzip DLL dependency test

### Changed
- Automate binary upload to GitHub Releases
- Bump version to 0.7.2 for release
- Consolidate build and test workflows into a single parallel pipeline
- Pysam search

### Fixed
- Add binaries resource path to PATH for DLL discovery

## [v0.7.1] - 2026-03-26
### Added
- Add pysam Linux verification

### Changed
- Bump version to 0.7.1
- Update Windows binary verification

## [v0.7.0] - 2026-03-25
### Added
- Add Windows verification system and fix samtools detection

### Changed
- Bump version to 0.7.0

### Fixed
- Fix test-windows workflow: use MSYS2 for sidecars
- Fix test-windows workflow: use conda for sidecars
- Fix test-windows workflow: use absolute paths for scoop shims
- Fix test-windows workflow: robust dependencies and use scoop for sidecars

## [v0.6.9] - 2026-03-25
### Fixed
- Revert to pip-only CI (pysam unavailable on Windows via conda)

## [v0.6.8] - 2026-03-25
### Fixed
- Remove htslib from conda install (bundled with pysam)

## [v0.6.7] - 2026-03-25
### Fixed
- Use conda for pysam on windows ci and clean resource PATH

## [v0.6.6] - 2026-03-25
### Fixed
- Include pysam in bundle and resolve race condition in job saving

## [v0.6.5] - 2026-03-25
### Fixed
- Use pysam for bgzip, add save retries, and support resource-path for DLLs

## [v0.6.4] - 2026-03-25
### Fixed
- Restore app initialization and improve sidecar detection

## [v0.6.3] - 2026-03-25
### Changed
- Trigger build only on releases

### Fixed
- Add CLI and environment logging

## [v0.6.2] - 2026-03-25
## [v0.6.1] - 2026-03-24
### Fixed
- Use robust subprocess.run without shell redirection for bgzip on Windows

## [v0.6.0] - 2026-03-24
### Added
- Implement automatic proxy detection and dynamic updates

## [v0.5.3] - 2026-03-23
### Added
- Samtools and bgzip paths
- VEP online and local

### Changed
- Dockerfile

## [v0.5.2] - 2026-02-26
### Added
- Add proxy support

## [v0.5.1] - 2026-02-26
### Added
- Enable searching reference by gene name

### Changed
- Update repository ownership URLs to lagosproject

## [v0.5.0] - 2026-02-24
### Fixed
- Enforce ubuntu-22.04 strictly for maximal linux backward/forward compatibility

## [v0.4.9] - 2026-02-24
### Fixed
- Bundle httpx properly by adding to root requirements.txt

## [v0.4.8] - 2026-02-24
### Fixed
- Change ubuntu-24.04 to ubuntu-latest in build matrix and bump version

## [v0.4.7] - 2026-02-24
### Fixed
- Add ubuntu-24.04 and ubuntu-22.04 to build matrix and append OS to output binary

## [v0.4.6] - 2026-02-24
### Fixed
- Downgrade ubuntu-latest runner to ubuntu-22.04 for GLIBC compatibility

## [v0.4.5] - 2026-02-24
### Fixed
- Bundle httpx for pyinstaller

## [v0.4.4] - 2026-02-24
### Fixed
- Install requirements individually to handle platform-specific failures

## [v0.4.3] - 2026-02-24
### Fixed
- Add site-packages to pathex, collect_submodules, and CI debug steps

## [v0.4.2] - 2026-02-24
### Fixed
- Add comprehensive hidden imports for uvicorn and dependencies in PyInstaller spec

## [v0.4.1] - 2026-02-24
### Added
- Delete unnecesary files

### Changed
- Trigger build on tag push

## [v0.4.0] - 2026-02-24
### Added
- Migrate to Ensembl/VEP batch HGVS annotator and fix RefSeq resolution

## [v0.3.0] - 2026-02-23
### Changed
- Rename MS Analyzer to PS Analyzer and bump version to 0.3.0
- Refactor to ps-analyzer

## [v0.2.1] - 2026-02-23
### Fixed
- Release workflow permissions and remove tracy bundling

## [v0.2.0] - 2026-02-23
### Added
- Import projects
- Tracy from path

### Changed
- Use .spec for builds and remove macos compilation
- Allow manual workflow dispatch
- Initial commit (Clean)
