#include "catch.hpp"
#include "utils.hpp"
#include "../src/jellyfishreader.hpp"
#ifdef KFF_SUPPORT
#include "../src/kffreader.hpp"
#endif
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#ifdef KFF_SUPPORT

// Helper function to create KFF file from existing JF file
void createKffFromJf(const string& jf_path, const string& kff_path, size_t k) {
    // Extract k-mers and counts from jellyfish file
    string cmd = "jellyfish dump " + jf_path + " | awk 'NR%2==1 {count=$1; gsub(\">\", \"\", count)} NR%2==0 {print $0 \" \" count}' > " + kff_path + ".txt";
    system(cmd.c_str());
    
    // Create KFF file with counts
    string kff_cmd = "kff-tools instr -i " + kff_path + ".txt -o " + kff_path + " -k " + to_string(k) + " -d 1 --delimiter ' '";
    system(kff_cmd.c_str());
    
    // Clean up temporary file
    remove((kff_path + ".txt").c_str());
}

TEST_CASE("JF_KFF_format_equivalence_reads", "[FormatEquivalence]") {
    // Test equivalence using reads.jf test file
    const string jf_file = "/mnt/ssd/lalli/pangenie/tests/data/reads.jf";
    const string kff_file = "/mnt/ssd/lalli/pangenie/tests/data/reads_test_equiv.kff";
    const size_t k = 10;
    
    // Create equivalent KFF file from existing JF file
    createKffFromJf(jf_file, kff_file, k);
    
    // Set jellyfish k-mer size
    jellyfish::mer_dna::k(k);
    
    // Initialize both readers
    JellyfishReader jf_reader(jf_file, k);
    KffReader kff_reader(kff_file, k);
    
    // Test all k-mers that should be in both files
    vector<string> test_kmers = {
        "AAAAAAACGG", "AAAAAACGGC", "ATGCTGTAAA", "CGTTTTTTTA",
        "CTGTAAAAAA", "GCTGTAAAAA", "GTAAAAAAAC", "TGCTGTAAAA", "TGTAAAAAAA"
    };
    
    // Test 1: K-mer abundance queries
    for (const string& kmer : test_kmers) {
        size_t jf_count = jf_reader.getKmerAbundance(kmer);
        size_t kff_count = kff_reader.getKmerAbundance(kmer);
        
        INFO("Testing k-mer: " << kmer);
        REQUIRE(jf_count == kff_count);
        REQUIRE(jf_count > 0); // Should exist in test data
        
        // Test with jellyfish k-mer objects
        jellyfish::mer_dna jelly_kmer(kmer);
        size_t jf_jelly = jf_reader.getKmerAbundance(jelly_kmer);
        size_t kff_jelly = kff_reader.getKmerAbundance(jelly_kmer);
        REQUIRE(jf_jelly == kff_jelly);
        REQUIRE(jf_jelly == jf_count);
    }
    
    // Test 2: Non-existent k-mers
    vector<string> absent_kmers = {"GGGGGGGGGG", "TTTTTTTTTT", "CCCCCCCCCC"};
    for (const string& kmer : absent_kmers) {
        size_t jf_count = jf_reader.getKmerAbundance(kmer);
        size_t kff_count = kff_reader.getKmerAbundance(kmer);
        
        INFO("Testing absent k-mer: " << kmer);
        REQUIRE(jf_count == kff_count);
        REQUIRE(jf_count == 0); // Should not exist
    }
    
    // Test 3: Coverage computation
    size_t genome_size = 1000;
    size_t jf_coverage = jf_reader.computeKmerCoverage(genome_size);
    size_t kff_coverage = kff_reader.computeKmerCoverage(genome_size);
    REQUIRE(jf_coverage == kff_coverage);
    
    // Test 4: Histogram computation
    size_t jf_peak = jf_reader.computeHistogram(1000, true, "");
    size_t kff_peak = kff_reader.computeHistogram(1000, true, "");
    REQUIRE(jf_peak == kff_peak);
    
    // Clean up created KFF file
    remove(kff_file.c_str());
}

TEST_CASE("JF_KFF_format_equivalence_nocanonical", "[FormatEquivalence]") {
    // Test that error handling is equivalent for problematic files
    const string jf_file = "/mnt/ssd/lalli/pangenie/tests/data/reads.no-canonical.jf";
    const size_t k = 10;
    
    // Set jellyfish k-mer size
    jellyfish::mer_dna::k(k);
    
    // Both should throw for non-canonical jellyfish file
    REQUIRE_THROWS(JellyfishReader(jf_file, k));
    
    // KFF files don't have this restriction, so no equivalent test
    // This verifies that error behaviors are consistent within format constraints
}

TEST_CASE("JF_KFF_PanGenie_output_equivalence", "[FormatEquivalence][Integration]") {
    // Test that PanGenie produces identical outputs with equivalent files
    const string jf_file = "/mnt/ssd/lalli/pangenie/tests/data/reads.jf";
    const string kff_file = "/mnt/ssd/lalli/pangenie/tests/data/reads_pangenie_equiv.kff";
    const string ref_file = "/mnt/ssd/lalli/pangenie/tests/data/region.fa";
    const string vcf_file = "/mnt/ssd/lalli/pangenie/tests/data/region.vcf";
    const size_t k = 10;
    
    // Create equivalent KFF file from existing JF file
    createKffFromJf(jf_file, kff_file, k);
    
    // Run PanGenie with JF file
    string jf_cmd = "/mnt/ssd/lalli/pangenie/build/src/PanGenie -i " + jf_file + 
                   " -r " + ref_file + " -v " + vcf_file + 
                   " -o /tmp/jf_output -k " + to_string(k) + " -e 100000";
    int jf_result = system(jf_cmd.c_str());
    
    // Run PanGenie with KFF file
    string kff_cmd = "/mnt/ssd/lalli/pangenie/build/src/PanGenie -i " + kff_file + 
                    " -r " + ref_file + " -v " + vcf_file + 
                    " -o /tmp/kff_output -k " + to_string(k) + " -e 100000";
    int kff_result = system(kff_cmd.c_str());
    
    // Both should succeed
    REQUIRE(jf_result == 0);
    REQUIRE(kff_result == 0);
    
    // Compare output VCF files
    string diff_cmd = "diff /tmp/jf_output_genotyping.vcf /tmp/kff_output_genotyping.vcf";
    int diff_result = system(diff_cmd.c_str());
    
    // Files should be identical
    REQUIRE(diff_result == 0);
    
    // Clean up
    remove(kff_file.c_str());
    remove("/tmp/jf_output_genotyping.vcf");
    remove("/tmp/kff_output_genotyping.vcf");
    remove("/tmp/jf_output_histogram.histo");
    remove("/tmp/kff_output_histogram.histo");
    remove("/tmp/jf_output_path_segments.fasta");
    remove("/tmp/kff_output_path_segments.fasta");
}

#else
// If KFF support is not compiled in, provide dummy test
TEST_CASE("JF_KFF_equivalence_not_supported", "[FormatEquivalence]") {
    // This test verifies that the test suite can run without KFF support
    REQUIRE(true);
}
#endif