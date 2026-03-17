# Changelog

All notable changes to ShigaPass will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.6.0] - 2026-03-17

### Added
- **Parallel processing support** with new `-j` option to process multiple samples simultaneously
- **Enhanced logging system** with individual sample log files and timestamps
- **Organized output structure** with dedicated `logs/` directory for sample-specific logs
- **Thread-safe file operations** using `flock` for concurrent CSV writing
- **Progress tracking** with clean sample-by-sample status updates
- **Performance optimization** guidelines in README for optimal thread/job configuration

### Changed
- **Input file format** now supports two-column format: `<Strain_ID> <path_to_contigs>` (backward compatible)
- **Log file organization** - main log stored in results directory root, individual sample logs in `logs/` subdirectory
- **Usage documentation** updated to include parallel processing options and performance guidelines
- **Output Files documentation** enhanced to describe new logging structure

### Improved
- **Stdout output** cleaned up - removed verbose analysis messages, keeping only essential progress indicators
- **Error handling** improved for parallel processing scenarios
- **Documentation** updated with parallel processing examples and best practices

### Technical Details
- Main log file: `ShigaPass_YYYYMMDD_HHMMSS.log` (stored in results directory)
- Individual sample logs: `logs/{strain_id}_ShigaPass.log` (timestamped entries with sample ID prefixes)
- Parallel job control with configurable limits (1-20 jobs, default: 1 for sequential processing)
- Thread allocation per BLAST operation optimized for parallel execution

## [1.5.0] - Previous Version

### Features
- Basic sequential sample processing
- Shigella serotype prediction using rfb, MLST, fliC, CRISPR, and ipaH markers
- Support for differentiating between Shigella, EIEC, and non-Shigella/EIEC
- Database initialization with `-u` option
- Optional intermediate file retention with `-k` option
- Comprehensive output with summary CSV files