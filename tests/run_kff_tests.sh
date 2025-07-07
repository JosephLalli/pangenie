#!/bin/bash

# Script to run KFF-related tests for PanGenie
# This script helps set up and run the KFF test suite

set -e

echo "=== PanGenie KFF Test Suite ==="

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    echo "Error: Please run this script from the tests directory"
    exit 1
fi

BUILD_DIR="../build"
DATA_DIR="data"

echo "1. Checking build directory..."
if [ ! -d "$BUILD_DIR" ]; then
    echo "Build directory not found. Please build PanGenie first:"
    echo "  mkdir ../build && cd ../build"
    echo "  cmake .. && make"
    exit 1
fi

echo "2. Checking for KFF support in build..."
if [ ! -f "$BUILD_DIR/tests/create_minimal_kff" ]; then
    echo "KFF support not found in build."
    echo "KFF tests will be skipped (this is expected if KFF C++ API is not installed)."
    KFF_AVAILABLE=false
else
    echo "KFF support detected!"
    KFF_AVAILABLE=true
fi

echo "3. Setting up test data..."
if [ "$KFF_AVAILABLE" = true ]; then
    echo "Creating test KFF files..."
    
    # Create test KFF file if it doesn't exist
    if [ ! -f "$DATA_DIR/test.kff" ]; then
        echo "  Creating test.kff..."
        $BUILD_DIR/tests/create_minimal_kff $DATA_DIR/test.kff
    else
        echo "  test.kff already exists"
    fi
    
    # Create additional test files for error testing
    if [ ! -f "$DATA_DIR/test-k15.kff" ]; then
        echo "  Creating test-k15.kff for error testing..."
        # Create a simple KFF file with k=15 for testing k-mer size mismatches
        cat > $DATA_DIR/temp_k15.txt << 'EOF'
ATGCTGTAAAAAAAC	1
TGCTGTAAAAAAACC	1
GCTGTAAAAAACCGG	1
EOF
        if command -v kff-tools >/dev/null 2>&1; then
            kff-tools instr -i $DATA_DIR/temp_k15.txt -o $DATA_DIR/test-k15.kff -k 15 2>/dev/null || true
            rm -f $DATA_DIR/temp_k15.txt
        else
            # If kff-tools not available, just create a dummy file
            echo "KFF" > $DATA_DIR/test-k15.kff
        fi
    fi
else
    echo "  Skipping KFF file creation (KFF support not available)"
fi

echo "4. Running tests..."
cd $BUILD_DIR/tests

if [ "$KFF_AVAILABLE" = true ]; then
    echo "Running all tests including KFF tests..."
    ./tests "[KFF]" -v
    echo ""
    echo "Running KFF-specific test cases..."
    ./tests "[KffReader]" -v
else
    echo "Running tests without KFF (KFF tests will be skipped)..."
    ./tests "~[KFF]" -v
fi

echo ""
echo "=== Test Summary ==="
echo "KFF Support: $([[ $KFF_AVAILABLE == true ]] && echo "Available" || echo "Not Available")"
echo "Test Status: Completed"

if [ "$KFF_AVAILABLE" = false ]; then
    echo ""
    echo "To enable KFF tests:"
    echo "1. Install KFF C++ API: https://github.com/Kmer-File-Format/kff-cpp-api"
    echo "2. Rebuild PanGenie with: cd ../build && cmake -DENABLE_KFF_SUPPORT=ON .. && make"
    echo "3. Run this script again"
fi