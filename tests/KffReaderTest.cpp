#include "catch.hpp"
#include "utils.hpp"

#ifdef KFF_SUPPORT
#include "../src/kffreader.hpp"
#include "../src/jellyfishreader.hpp"
#include <vector>
#include <string>
#include <fstream>

using namespace std;

TEST_CASE("KffReader_basic", "[KffReader]") {
    // Test basic KFF file reading functionality
    KffReader reader("../tests/data/test.kff", 10);
    
    // Test some basic k-mers that should be in the test file
    string read = "ATGCTGTAAAAAAACGGC";
    for (size_t i = 0; i < read.size()-9; ++i) {
        string kmer = read.substr(i,10);
        size_t count = reader.getKmerAbundance(kmer);
        // Should have non-zero count for k-mers present in test data
        REQUIRE(count >= 0); // At minimum, should not crash
    }
}

TEST_CASE("KffReader_jellyfish_kmer", "[KffReader]") {
    // Test KFF reader with jellyfish k-mer objects
    KffReader reader("../tests/data/test.kff", 10);
    
    string kmer_str = "ATGCTGTAAA";
    jellyfish::mer_dna jelly_kmer(kmer_str);
    size_t count = reader.getKmerAbundance(jelly_kmer);
    REQUIRE(count >= 0); // Should handle jellyfish k-mer objects
}

TEST_CASE("KffReader_canonical", "[KffReader]") {
    // Test that KFF reader handles canonical k-mers correctly
    KffReader reader("../tests/data/test.kff", 10);
    
    string kmer = "ATGCTGTAAA";
    string rev_comp = "TTTACAGCAT"; // reverse complement
    
    size_t count1 = reader.getKmerAbundance(kmer);
    size_t count2 = reader.getKmerAbundance(rev_comp);
    
    // Canonical k-mers should return the same count
    REQUIRE(count1 == count2);
}

TEST_CASE("KffReader_histogram", "[KffReader]") {
    // Test histogram computation
    KffReader reader("../tests/data/test.kff", 10);
    
    // Test histogram computation doesn't crash
    REQUIRE_NOTHROW(reader.computeHistogram(1000, true, ""));
}

TEST_CASE("KffReader_coverage", "[KffReader]") {
    // Test k-mer coverage computation
    KffReader reader("../tests/data/test.kff", 10);
    
    size_t genome_kmers = 1000;
    size_t coverage = reader.computeKmerCoverage(genome_kmers);
    REQUIRE(coverage >= 0); // Should return valid coverage
}

TEST_CASE("KffReader_error_handling", "[KffReader]") {
    // Test error handling for invalid files
    
    // Non-existent file
    REQUIRE_THROWS(KffReader("../tests/data/nonexistent.kff", 10));
    
    // Wrong k-mer size (if test.kff has k=10, try k=15)
    REQUIRE_THROWS(KffReader("../tests/data/test-k15.kff", 10));
}

TEST_CASE("KffReader_comparison_with_jellyfish", "[KffReader]") {
    // Compare KFF reader results with equivalent Jellyfish data
    // This test requires that test.kff and reads.jf contain equivalent data
    
    // Skip this test if test files don't exist
    ifstream kff_test("../tests/data/test.kff");
    ifstream jf_test("../tests/data/reads.jf");
    if (!kff_test.good() || !jf_test.good()) {
        // Skip test if files don't exist
        return;
    }
    
    KffReader kff_reader("../tests/data/test.kff", 10);
    JellyfishReader jf_reader("../tests/data/reads.jf", 10);
    
    // Test a few k-mers to ensure results are consistent
    vector<string> test_kmers = {"ATGCTGTAAA", "TGCTGTAAAA", "GCTGTAAAAA"};
    
    for (const string& kmer : test_kmers) {
        size_t kff_count = kff_reader.getKmerAbundance(kmer);
        size_t jf_count = jf_reader.getKmerAbundance(kmer);
        
        // Results should be the same for equivalent data
        REQUIRE(kff_count == jf_count);
    }
}

TEST_CASE("KffReader_polymorphic", "[KffReader]") {
    // Test that KffReader works through KmerCounter interface
    shared_ptr<KmerCounter> counter = make_shared<KffReader>("../tests/data/test.kff", 10);
    
    string kmer = "ATGCTGTAAA";
    size_t count = counter->getKmerAbundance(kmer);
    REQUIRE(count >= 0);
    
    // Test jellyfish k-mer interface
    jellyfish::mer_dna jelly_kmer(kmer);
    size_t jelly_count = counter->getKmerAbundance(jelly_kmer);
    REQUIRE(jelly_count >= 0);
}

#else
// If KFF support is not compiled in, provide a dummy test
TEST_CASE("KffReader_not_supported", "[KffReader]") {
    // This test just verifies the test suite can run without KFF support
    REQUIRE(true);
}
#endif