# KFF Support in PanGenie

This document describes how to enable and use KFF (K-mer File Format) support in PanGenie.

## Overview

PanGenie now supports KFF files as an alternative to Jellyfish (.jf) files for pre-computed k-mer counts. KFF is a standardized, space-efficient format for storing k-mers and their associated data.

## Benefits of KFF Format

- **Standardized**: Universal format supported by multiple tools
- **Space-efficient**: ~3x more compact than native formats
- **Interoperability**: Works with KMC, DSK, kmtricks, and other tools
- **Future-proof**: Growing ecosystem support

## Installation

### Prerequisites

1. **KFF C++ API**: Required for KFF support
```bash
git clone https://github.com/Kmer-File-Format/kff-cpp-api
cd kff-cpp-api
mkdir build && cd build
cmake .. && make && sudo make install
```

2. **Standard PanGenie dependencies**: jellyfish, cereal, zlib, cmake, gcc

### Building with KFF Support

```bash
# Clone and setup PanGenie
git clone https://github.com/eblerjana/pangenie.git
cd pangenie
conda env create -f environment.yml
conda activate pangenie

# Build with KFF support enabled
mkdir build && cd build
cmake -DENABLE_KFF_SUPPORT=ON ..
make

# Verify KFF support was enabled
# Look for "Found KFF C++ API" message during cmake
```

### Building without KFF Support (default)

```bash
# Standard build (KFF support disabled by default)
mkdir build && cd build
cmake ..
make
```

## Usage

### Using KFF Files

When KFF support is enabled, PanGenie automatically detects `.kff` files:

```bash
# Using KFF file with preprocessing workflow
PanGenie-index -r reference.fa -v variants.vcf -o preprocessing
PanGenie -f preprocessing -i reads.kff -o results

# Using KFF file in single-step mode
PanGenie -i reads.kff -r reference.fa -v variants.vcf -s sample_name
```

### Creating KFF Files

You can create KFF files using various tools:

#### Using KMC + kff-tools
```bash
# Count k-mers with KMC
kmc -k31 -m8 reads.fa reads_kmc ./temp/

# Convert to KFF format
kmc_dump reads_kmc reads.txt
kff-tools instr -i reads.txt -o reads.kff -k 31
```

#### Using DSK (with KFF output)
```bash
# DSK can directly output KFF format
dsk -file reads.fa -kmer-size 31 -out reads.kff -out-fmt kff
```

## CMake Configuration Options

```bash
# Enable KFF support (fails if KFF API not found)
cmake -DENABLE_KFF_SUPPORT=ON ..

# Disable KFF support explicitly
cmake -DENABLE_KFF_SUPPORT=OFF ..

# Default behavior (KFF disabled)
cmake ..
```

## Testing

### Running KFF Tests

```bash
# Run all tests including KFF (if enabled)
cd build/tests
make check

# Run only KFF-specific tests
./tests "[KFF]" -v

# Automated test setup and execution
cd tests
./run_kff_tests.sh
```

### Test Coverage

KFF tests include:
- Basic file reading and k-mer lookup
- Canonical k-mer handling
- Histogram and coverage computation  
- Error handling and edge cases
- Integration with commands.cpp
- Compatibility with existing Jellyfish workflows

## Troubleshooting

### CMake Issues

**Error: "KFF C++ API not found"**
```bash
# Solution: Install KFF C++ API or disable support
cmake -DENABLE_KFF_SUPPORT=OFF ..
```

**Error: "kff_io.hpp not found"**
```bash
# Ensure KFF headers are installed in standard locations
sudo ldconfig  # Refresh library cache
# Or specify custom paths:
cmake -DENABLE_KFF_SUPPORT=ON -DKFF_INCLUDE_DIR=/custom/path/include ..
```

### Runtime Issues

**Error: "Unknown file format"**
- Verify file is valid KFF format
- Check k-mer size matches between file and PanGenie parameters
- Ensure file is not corrupted

**Performance**: KFF files use in-memory caching for fast lookups. Large files may require significant RAM.

## Backward Compatibility

- Existing Jellyfish (.jf) workflows continue to work unchanged
- FASTA/FASTQ input processing remains the same
- All PanGenie features work with both file formats
- No changes required for existing scripts/pipelines

## File Format Comparison

| Feature | Jellyfish (.jf) | KFF (.kff) |
|---------|----------------|------------|
| Space efficiency | Baseline | ~3x smaller |
| Tool support | Jellyfish only | Multiple tools |
| Standardization | Tool-specific | Universal |
| PanGenie support | Native | Optional |
| Random access | Yes | Sequential |
| Creation speed | Fast | Moderate |

## Contributing

When contributing KFF-related code:
- Always use `#ifdef KFF_SUPPORT` guards
- Ensure tests work with and without KFF support
- Update documentation for new features
- Test both enabled and disabled build configurations