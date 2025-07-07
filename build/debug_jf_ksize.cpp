#include <iostream>
#include <fstream>
#include <stdexcept>
#include <jellyfish/file_header.hpp>
#include <jellyfish/jellyfish.hpp>

using namespace std;

int main() {
    try {
        // Open the jellyfish file
        string filename = "../tests/data/reads.jf";
        ifstream ifs(filename, ios::in | ios::binary);
        
        if (!ifs.good()) {
            cerr << "ERROR: Cannot open file: " << filename << endl;
            return 1;
        }
        
        // Read the file header
        jellyfish::file_header header(ifs);
        
        if (!ifs.good()) {
            cerr << "ERROR: Failed to parse header of file: " << filename << endl;
            return 1;
        }
        
        // Print metadata
        cout << "=== Jellyfish File Metadata ===" << endl;
        cout << "File: " << filename << endl;
        cout << "Key length (bits): " << header.key_len() << endl;
        cout << "K-mer size: " << (header.key_len() / 2) << endl;
        cout << "Counter length: " << header.counter_len() << endl;
        cout << "Canonical: " << (header.canonical() ? "Yes" : "No") << endl;
        cout << "Format: " << header.format() << endl;
        cout << "Size: " << header.size() << endl;
        cout << "Offset: " << header.offset() << endl;
        cout << "Matrix rows: " << header.matrix().r() << endl;
        cout << "Matrix cols: " << header.matrix().c() << endl;
        cout << "================================" << endl;
        
        // Check if format is supported
        if (header.format() == binary_dumper::format) {
            cout << "Format is binary dumper (supported)" << endl;
        } else {
            cout << "Format is NOT binary dumper (unsupported)" << endl;
        }
        
        ifs.close();
        
    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}