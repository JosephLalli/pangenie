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