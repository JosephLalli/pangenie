#include <iostream>
#include "../src/jellyfishreader.hpp"
#include "../src/kffreader.hpp"
#include <jellyfish/mer_dna.hpp>

int main() {
    try {
        jellyfish::mer_dna::k(10);
        
        JellyfishReader jf_reader("../tests/data/reads.jf", 10);
        KffReader kff_reader("../tests/data/reads_equivalent.kff", 10);
        
        std::string test_kmer = "ATGCTGTAAA";
        std::cout << "Testing kmer: " << test_kmer << std::endl;
        std::cout << "JF count: " << jf_reader.getKmerAbundance(test_kmer) << std::endl;
        std::cout << "KFF count: " << kff_reader.getKmerAbundance(test_kmer) << std::endl;
        
        // Test coverage computation
        std::cout << "JF coverage: " << jf_reader.computeKmerCoverage(1000) << std::endl;
        std::cout << "KFF coverage: " << kff_reader.computeKmerCoverage(1000) << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}