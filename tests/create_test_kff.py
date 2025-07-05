#!/usr/bin/env python3
"""
Script to create test KFF files for PanGenie testing.
This script creates KFF files from the existing test data to ensure 
compatibility with the existing test suite.
"""

import subprocess
import sys
import os

def create_kff_from_fasta(fasta_file, kff_file, k=10):
    """
    Create a KFF file from a FASTA file using kff-tools.
    
    Args:
        fasta_file: Path to input FASTA file
        kff_file: Path to output KFF file
        k: K-mer size (default: 10)
    """
    try:
        # Use kmc to count k-mers first, then convert to KFF
        temp_prefix = kff_file.replace('.kff', '')
        
        # Step 1: Create k-mer database with KMC
        cmd_kmc = [
            'kmc', 
            '-k{}'.format(k), 
            '-m4',  # 4GB memory limit
            '-sm',  # small k-mer mode
            '-ci1', # minimum count = 1
            '-cs100000', # maximum count = 100000
            fasta_file,
            temp_prefix,
            '.'  # temporary directory
        ]
        
        print(f"Creating k-mer database: {' '.join(cmd_kmc)}")
        result = subprocess.run(cmd_kmc, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"KMC failed: {result.stderr}")
            return False
        
        # Step 2: Convert KMC database to KFF using kff-tools
        cmd_convert = [
            'kmc_dump',
            temp_prefix,
            temp_prefix + '.txt'
        ]
        
        print(f"Dumping k-mers: {' '.join(cmd_convert)}")
        result = subprocess.run(cmd_convert, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"KMC dump failed: {result.stderr}")
            return False
        
        # Step 3: Convert text to KFF using kff-tools
        cmd_kff = [
            'kff-tools',
            'instr',
            '-i', temp_prefix + '.txt',
            '-o', kff_file,
            '-k', str(k)
        ]
        
        print(f"Creating KFF file: {' '.join(cmd_kff)}")
        result = subprocess.run(cmd_kff, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"KFF creation failed: {result.stderr}")
            return False
        
        # Clean up temporary files
        for ext in ['.kmc_pre', '.kmc_suf', '.txt']:
            temp_file = temp_prefix + ext
            if os.path.exists(temp_file):
                os.remove(temp_file)
        
        print(f"Successfully created {kff_file}")
        return True
        
    except FileNotFoundError as e:
        print(f"Required tool not found: {e}")
        print("Please install KMC and kff-tools")
        return False
    except Exception as e:
        print(f"Error creating KFF file: {e}")
        return False

def create_simple_kff(kff_file, k=10):
    """
    Create a simple KFF file with known k-mers for testing.
    This creates a KFF file directly using the KFF API if available.
    """
    try:
        # Create a simple text file with k-mers and counts
        temp_txt = kff_file.replace('.kff', '.txt')
        
        # Known k-mers from the test data
        test_kmers = [
            ("ATGCTGTAAA", 1),
            ("TGCTGTAAAA", 1), 
            ("GCTGTAAAAA", 1),
            ("CTGTAAAAAA", 1),
            ("TGTAAAAAAA", 1),
            ("GTAAAAACGG", 1),
            ("TAAAAACGGC", 1)
        ]
        
        with open(temp_txt, 'w') as f:
            for kmer, count in test_kmers:
                f.write(f"{kmer}\t{count}\n")
        
        # Convert to KFF using kff-tools
        cmd_kff = [
            'kff-tools',
            'instr',
            '-i', temp_txt,
            '-o', kff_file,
            '-k', str(k)
        ]
        
        print(f"Creating simple KFF file: {' '.join(cmd_kff)}")
        result = subprocess.run(cmd_kff, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"KFF creation failed: {result.stderr}")
            return False
        
        # Clean up
        if os.path.exists(temp_txt):
            os.remove(temp_txt)
        
        print(f"Successfully created simple {kff_file}")
        return True
        
    except Exception as e:
        print(f"Error creating simple KFF file: {e}")
        return False

def main():
    """Create test KFF files for the test suite."""
    
    test_data_dir = os.path.join(os.path.dirname(__file__), 'data')
    
    # Check if test data directory exists
    if not os.path.exists(test_data_dir):
        print(f"Test data directory not found: {test_data_dir}")
        return 1
    
    # Files to create
    kff_files = [
        ('test.kff', 10),        # Basic test file
        ('test-k15.kff', 15),    # Different k-mer size for error testing
    ]
    
    success_count = 0
    
    for kff_filename, k in kff_files:
        kff_path = os.path.join(test_data_dir, kff_filename)
        
        print(f"\nCreating {kff_filename} with k={k}...")
        
        # Try to create from existing FASTA file first
        fasta_path = os.path.join(test_data_dir, 'reads.fa')
        if os.path.exists(fasta_path) and k == 10:
            if create_kff_from_fasta(fasta_path, kff_path, k):
                success_count += 1
                continue
        
        # Fall back to creating simple KFF file
        if create_simple_kff(kff_path, k):
            success_count += 1
        else:
            print(f"Failed to create {kff_filename}")
    
    print(f"\nCreated {success_count}/{len(kff_files)} KFF test files")
    
    if success_count == 0:
        print("\nNo KFF files created. Tests will be skipped.")
        print("To enable KFF tests, install KMC and kff-tools:")
        print("  - KMC: https://github.com/refresh-bio/KMC")
        print("  - kff-tools: https://github.com/Kmer-File-Format/kff-tools")
    
    return 0 if success_count > 0 else 1

if __name__ == '__main__':
    sys.exit(main())