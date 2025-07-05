# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PanGenie is a short-read genotyper for genetic variants (SNPs, indels, SVs) represented in pangenome graphs. It uses k-mer counting and known haplotype sequences to compute genotypes. The project is written in C++ and uses CMake for building.

## Build System & Commands

### Building the Project

```bash
# Using conda (recommended)
conda env create -f environment.yml
conda activate pangenie
mkdir build && cd build

# Basic build (without KFF support)
cmake ..
make

# Build with KFF support (requires KFF C++ API)
cmake -DENABLE_KFF_SUPPORT=ON ..
make

# Using Singularity (alternative)
sudo singularity build pangenie.sif container/pangenie.def
```

### Running Tests

```bash
# Build tests
make check

# Run tests manually
cd build
ctest
```

### Demo Usage

```bash
# Quick demo (from demo/ directory)
PanGenie-index -r test-reference.fa -v test-variants.vcf -o preprocessing -e 100000
PanGenie -f preprocessing -i test-reads.fa -o test -e 100000
```

## Dependencies

- **jellyfish** (v2.2.10+): K-mer counting library
- **cereal**: C++ serialization library
- **zlib**: Compression library
- **cmake** (v3.1+): Build system
- **gcc** (v11.2+): C++20 compiler
- **kff-cpp-api** (optional): KFF file format support

## Architecture

### Core Executables

1. **PanGenie-index**: Pre-processing step that creates k-mer indices from reference and VCF files
2. **PanGenie**: Main genotyping engine using k-mer counting and HMM algorithms
3. **PanGenie-vcf**: Converts serialized genotyping results to VCF format
4. **PanGenie-sampling**: Panel sampling functionality to reduce computational complexity
5. **Analyze-UK**: Analysis tool for unique k-mers debugging

### Key Libraries & Components

- **PanGenieLib**: Shared library containing core functionality
- **Graph**: Pangenome graph representation and operations
- **HMM**: Hidden Markov Model for genotyping and phasing
- **KmerPath**: K-mer path handling and unique k-mer computation
- **EmissionProbabilityComputer**: Probability calculations for genotyping
- **VariantReader**: VCF file parsing and variant processing

### Workflow Architecture

1. **Index Creation**: Process reference genome and VCF to create k-mer indices
2. **K-mer Counting**: Count k-mers in sequencing reads using jellyfish
3. **Genotyping**: Use Forward-Backward algorithm for genotype inference
4. **Phasing**: Optional Viterbi algorithm for haplotype phasing
5. **Output**: Generate VCF files with genotype calls

## File Structure

- `src/`: C++ source code and main executables
- `tests/`: Unit tests with catch.hpp framework
- `demo/`: Small test datasets for verification
- `scripts/`: Python analysis and utility scripts
- `pipelines/`: Snakemake workflows for callset processing
- `container/`: Singularity container definition

## Development Notes

### Testing

Tests use the catch.hpp framework. Test data is located in `tests/data/` directory. Each core component has corresponding test files (e.g., `GraphTest.cpp`, `HMMTest.cpp`).

### Memory Management

PanGenie uses serialization (cereal) for efficient memory management and data persistence. Index files are created once and reused across multiple samples.

### Threading

The tool supports multi-threading for both k-mer counting (`-j` parameter) and genotyping (`-t` parameter). Different chromosomes can be processed in parallel.

### Input Requirements

- **VCF files**: Must be multi-sample, fully-phased, non-overlapping variants, sequence-resolved
- **Reference**: Standard FASTA format (uncompressed)
- **Reads**: FASTA/FASTQ format (uncompressed), pre-computed jellyfish database (.jf), or KFF files (.kff)

## Code Patterns

### Error Handling

The codebase uses exceptions for error handling. Check existing patterns in `commandlineparser.cpp` and `variantreader.cpp`.

### Memory Allocation

Uses modern C++ patterns with smart pointers and RAII. Large data structures are serialized to disk for memory efficiency.

### K-mer Handling

K-mer operations are optimized using jellyfish integration. Default k-mer size is 31, configurable via `-k` parameter.

### KFF File Support

PanGenie supports KFF (K-mer File Format) files when explicitly enabled:
- **Optional feature**: Must be enabled with `-DENABLE_KFF_SUPPORT=ON` during cmake configuration
- **Automatic detection**: Files with `.kff` extension are automatically recognized when enabled
- **Conditional compilation**: KFF support is enabled via `KFF_SUPPORT` preprocessor flag
- **Backward compatibility**: Existing jellyfish `.jf` file support remains unchanged
- **Installation**: KFF C++ API must be installed and findable by CMake

#### Enabling KFF Support

1. Install KFF C++ API:
```bash
git clone https://github.com/Kmer-File-Format/kff-cpp-api
cd kff-cpp-api
mkdir build && cd build
cmake .. && make && sudo make install
```

2. Build PanGenie with KFF support:
```bash
cmake -DENABLE_KFF_SUPPORT=ON ..
make
```

3. Verify KFF support:
```bash
# Should show KFF-related messages during cmake configuration
# Tests will include KFF functionality when enabled
```

## Current Status

### ✅ Working
- **PanGenie builds successfully with KFF support enabled**
- **Main executables**: `PanGenie` and `PanGenie-index` compiled with KFF symbols
- **CMake configuration**: Properly detects KFF C++ API at `/home/lalli/usr/local`
- **Environment**: CMAKE_PREFIX_PATH correctly set to `/home/lalli/usr/local`

### ⚠️ Known Issues
- **Test compatibility**: Some test files may have compatibility issues with newer KFF API versions
- **Symbolic links**: Fixed broken KFF header links in `/home/lalli/usr/local/include/`

## Lessons Learned: CMake & Build Environment

### Local Installation Paths
- **Critical**: KFF C++ API is installed in `/home/lalli/usr/local`, NOT `/usr/local`
- **Environment**: CMAKE_PREFIX_PATH correctly set to `/home/lalli/usr/local`
- **CMake configuration**: Updated to use `/home/lalli/usr/local` as default install prefix

### Symbolic Link Management
- **Problem**: Broken symbolic links in `/home/lalli/usr/local/include/` caused header detection failures
- **Solution**: Fixed links to point to actual header files in build directory
- **Verification**: Use `ls -la` to check link validity before troubleshooting CMake

### CMake Debugging
- **Tool**: Use `cmake --debug-output` to trace exact search paths
- **Pattern**: CMake finds libraries but fails on headers due to broken links
- **Verification**: Check both `find_path` and `find_library` results separately

### Safety Practices
- **File operations**: Always use absolute paths with `rm` commands
- **Verification**: Confirm current directory before running destructive commands
- **Never reinstall**: Don't attempt to reinstall KFF C++ API without explicit permission

### Build Configuration
- **Threading**: Use `-j24` for parallel compilation
- **KFF Detection**: CMake lines 40-64 handle KFF detection logic
- **Verification**: Check binary symbols with `strings` command to confirm KFF integration

## Troubleshooting KFF Integration

1. **Check environment**: Verify CMAKE_PREFIX_PATH includes `/home/lalli/usr/local`
2. **Verify installation**: Confirm KFF library and headers exist in local install directory
3. **Check symbolic links**: Ensure header links point to valid files
4. **Debug CMake**: Use `--debug-output` to trace search paths
5. **Verify integration**: Check compiled binary contains KFF symbols

## KFF Test Data Creation

### Successfully Created KFF Test Files

**Completed**: Created `reads.kff` as KFF equivalent of existing `reads.jf` test data

#### Conversion Process Used:
```bash
# 1. Extract k-mers from existing jellyfish file
jellyfish dump reads.jf > reads_canonical_kmers.txt

# 2. Convert format (jellyfish dump uses >count\nkmer format)
awk 'NR%2==0' reads_canonical_kmers.txt > reads_kmer_only.txt  

# 3. Create KFF file using kff-tools
kff-tools instr -i reads_kmer_only.txt -o reads.kff -k 10
```

#### Validation Results:
- ✅ **File created**: `reads.kff` (209 bytes)
- ✅ **K-mer count**: 9 k-mers (matches expected from 18-line jellyfish dump)
- ✅ **Format verified**: `kff-tools outstr` confirms readable KFF format
- ✅ **Source parameters**: k=10, canonical=yes, from `reads.fa`

#### Next Steps Available:
- Test PanGenie functionality with the new KFF file
- Update test suite to include KFF file testing

#### Tools Required:
- `kff-tools` (available at `/home/lalli/usr/local/bin/kff-tools`)
- `jellyfish` for original data extraction

## Current Working Directory

Project root: `/mnt/ssd/lalli/pangenie`
Build directory: `/mnt/ssd/lalli/pangenie/build`