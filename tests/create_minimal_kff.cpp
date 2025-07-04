#ifdef KFF_SUPPORT
#include "kff_io.hpp"
#include <iostream>
#include <vector>
#include <string>

/**
 * Creates a minimal KFF file for testing purposes.
 * This creates a valid KFF file with known k-mers that can be used in tests.
 */
int create_test_kff(const std::string& filename, size_t k = 10) {
    try {
        // Create KFF file writer
        Kff_file kff_file(filename, "w");
        
        // Set global variables
        kff_file.write_var("k", (uint64_t)k);
        kff_file.write_var("data_size", (uint64_t)8); // 8 bytes for count
        kff_file.write_var("max", (uint64_t)255);
        
        // Define test k-mers and their counts
        std::vector<std::pair<std::string, uint64_t>> test_data = {
            {"ATGCTGTAAA", 1},
            {"TGCTGTAAAA", 1}, 
            {"GCTGTAAAAA", 1},
            {"CTGTAAAAAA", 1},
            {"TGTAAAAAAA", 1},
            {"GTAAAAACGG", 1},
            {"TAAAAACGGC", 1},
            {"AAAAACGGCC", 1}
        };
        
        // Create a raw section
        Section_Raw section(&kff_file);
        
        // Write k-mers to the section
        for (const auto& pair : test_data) {
            const std::string& kmer_str = pair.first;
            uint64_t count = pair.second;
            
            if (kmer_str.length() != k) {
                std::cerr << "Skipping k-mer " << kmer_str << " (wrong length)" << std::endl;
                continue;
            }
            
            // Convert string k-mer to binary
            uint8_t* kmer_binary = new uint8_t[(k + 3) / 4];
            memset(kmer_binary, 0, (k + 3) / 4);
            
            for (size_t i = 0; i < k; i++) {
                uint8_t nucleotide = 0;
                switch (kmer_str[i]) {
                    case 'A': case 'a': nucleotide = 0; break;
                    case 'C': case 'c': nucleotide = 1; break;
                    case 'G': case 'g': nucleotide = 2; break;
                    case 'T': case 't': nucleotide = 3; break;
                    default: nucleotide = 0; break;
                }
                
                uint8_t byte_idx = i / 4;
                uint8_t bit_offset = (i % 4) * 2;
                kmer_binary[byte_idx] |= (nucleotide << bit_offset);
            }
            
            // Write k-mer and count
            uint8_t* data = reinterpret_cast<uint8_t*>(&count);
            section.write_kmer(kmer_binary, data);
            
            delete[] kmer_binary;
        }
        
        section.close();
        kff_file.close();
        
        std::cout << "Created test KFF file: " << filename << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error creating KFF file: " << e.what() << std::endl;
        return 1;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <output.kff>" << std::endl;
        return 1;
    }
    
    return create_test_kff(argv[1]);
}

#else
#include <iostream>
int main() {
    std::cout << "KFF support not compiled in" << std::endl;
    return 1;
}
#endif