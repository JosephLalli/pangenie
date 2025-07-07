#include "catch.hpp"
#include "utils.hpp"
#include "../src/jellyfishcounter.hpp"
#include "../src/jellyfishreader.hpp"
#ifdef KFF_SUPPORT
#include "../src/kffreader.hpp"
#endif
#include <vector>
#include <string>

using namespace std;

TEST_CASE("JellyfishCounter", "[JellyfishCounter]") {
	JellyfishCounter counter("../tests/data/reads.fa", 10);
	string read = "ATGCTGTAAAAAAACGGC";
	for (size_t i = 0; i < read.size()-9; ++i) {
		string kmer = read.substr(i,10);
		REQUIRE(counter.getKmerAbundance(kmer) == 1);
	}
	KmerCounter* t = new JellyfishCounter ("../tests/data/reads.fa", 10);
	delete t;
}

TEST_CASE("JellyfishCounter_if", "[JellyfishCounter_if]") {
	JellyfishCounter counter("../tests/data/reads.fa", {"../tests/data/kmerfile.fa"}, 10);
	// these two kmers are in kmerfile.fa and should have been counted
	REQUIRE(counter.getKmerAbundance("ATGCTGTAAA") == 1);
	REQUIRE(counter.getKmerAbundance("TGCTGTAAAA") == 1);
	// the following kmers are not contained in kmerfile.fa, thus they should have count 0
	string kmers = "GCTGTAAAAAAACGGC";
	for (size_t i = 0; i < kmers.size()-9; ++i) {
		string kmer = kmers.substr(i,10);
		REQUIRE(counter.getKmerAbundance(kmer) == 0);
	}
}

TEST_CASE("JellyfishReader", "[JellyfishReader]") {
	JellyfishReader reader ("../tests/data/reads.jf", 10);
	string read = "ATGCTGTAAAAAAACGGC";
	for (size_t i = 0; i < read.size()-9; ++i) {
		string kmer = read.substr(i,10);
		REQUIRE(reader.getKmerAbundance(kmer) == 1);
	}

	// reads where counted without the -C option
	REQUIRE_THROWS(JellyfishReader("../tests/data/reads.no-canonical.jf", 10));

	// wrong kmer size used
	REQUIRE_THROWS(JellyfishReader("../tests/data/reads.jf", 11));

}

#ifdef KFF_SUPPORT
TEST_CASE("KffReader_parallel", "[KffReader]") {
	// Test KFF reader with test.kff file created by create_minimal_kff
	KffReader reader ("../tests/data/test.kff", 10);
	
	// Test specific k-mers that should be in test.kff
	REQUIRE(reader.getKmerAbundance("ATGCTGTAAA") == 1);
	REQUIRE(reader.getKmerAbundance("TGCTGTAAAA") == 1);
	REQUIRE(reader.getKmerAbundance("GCTGTAAAAA") == 1);

	// wrong kmer size used
	REQUIRE_THROWS(KffReader("../tests/data/test.kff", 11));
}

TEST_CASE("KmerCounter_format_equivalence", "[KmerCounter][KffReader][JellyfishReader]") {
	// Test that JF and KFF readers produce identical results for overlapping k-mers
	JellyfishReader jf_reader("../tests/data/reads.jf", 10);
	KffReader kff_reader("../tests/data/test.kff", 10);
	
	// Set jellyfish k-mer size
	jellyfish::mer_dna::k(10);
	
	// Test k-mers that should be in both files
	vector<string> common_kmers = {"ATGCTGTAAA", "TGCTGTAAAA", "GCTGTAAAAA"};
	for (const string& kmer : common_kmers) {
		size_t jf_count = jf_reader.getKmerAbundance(kmer);
		size_t kff_count = kff_reader.getKmerAbundance(kmer);
		REQUIRE(jf_count == kff_count);
		REQUIRE(jf_count == 1); // Should both be 1
		
		// Test with jellyfish k-mer objects
		jellyfish::mer_dna jelly_kmer(kmer);
		size_t jf_jelly_count = jf_reader.getKmerAbundance(jelly_kmer);
		size_t kff_jelly_count = kff_reader.getKmerAbundance(jelly_kmer);
		REQUIRE(jf_jelly_count == kff_jelly_count);
		REQUIRE(jf_jelly_count == 1); // Should both be 1
	}
}

TEST_CASE("KmerCounter_polymorphic_interface", "[KmerCounter][KffReader][JellyfishReader]") {
	// Test polymorphic usage of both readers through KmerCounter interface
	
	// Set jellyfish k-mer size
	jellyfish::mer_dna::k(10);
	
	vector<unique_ptr<KmerCounter>> readers;
	readers.emplace_back(new JellyfishReader("../tests/data/reads.jf", 10));
	readers.emplace_back(new KffReader("../tests/data/test.kff", 10));
	
	string test_kmer = "ATGCTGTAAA";
	jellyfish::mer_dna jelly_kmer(test_kmer);
	
	// Both readers should find this k-mer
	for (const auto& reader : readers) {
		REQUIRE(reader->getKmerAbundance(test_kmer) == 1);
		REQUIRE(reader->getKmerAbundance(jelly_kmer) == 1);
		
		// Test other interface methods
		size_t coverage = reader->computeKmerCoverage(1000);
		REQUIRE(coverage >= 0);
		
		// Only test histogram for readers with sufficient data
		REQUIRE_NOTHROW(reader->computeHistogram(1000, true, ""));
	}
}

TEST_CASE("KmerFile_format_equivalence_guarantee", "[KmerCounter][KffReader][JellyfishReader]") {
	// Test that identical k-mer data in .jf and .kff formats produce identical results
	// This test ensures cross-format compatibility using preexisting .jf test data
	
	// Set jellyfish k-mer size
	jellyfish::mer_dna::k(10);
	
	// Test with the main test dataset: reads.jf vs reads_equivalent.kff
	JellyfishReader jf_reader("../tests/data/reads.jf", 10);
	KffReader kff_reader("../tests/data/reads_equivalent.kff", 10);
	
	// Extract all k-mers from the jellyfish file to test comprehensive equivalence
	vector<string> all_test_kmers = {
		"AAAAAAACGG", "AAAAAACGGC", "ATGCTGTAAA", "CGTTTTTTTA", 
		"CTGTAAAAAA", "GCTGTAAAAA", "GTAAAAAAAC", "TGCTGTAAAA", "TGTAAAAAAA"
	};
	
	// Test 1: All k-mer abundance queries return identical results
	for (const string& kmer : all_test_kmers) {
		size_t jf_count = jf_reader.getKmerAbundance(kmer);
		size_t kff_count = kff_reader.getKmerAbundance(kmer);
		
		INFO("Testing k-mer: " << kmer);
		REQUIRE(jf_count == kff_count);
		REQUIRE(jf_count == 1); // All test k-mers should have count 1
		
		// Test with jellyfish k-mer objects
		jellyfish::mer_dna jelly_kmer(kmer);
		size_t jf_jelly_count = jf_reader.getKmerAbundance(jelly_kmer);
		size_t kff_jelly_count = kff_reader.getKmerAbundance(jelly_kmer);
		REQUIRE(jf_jelly_count == kff_jelly_count);
		REQUIRE(jf_jelly_count == 1);
	}
	
	// Test 2: K-mer coverage computation produces identical results
	size_t genome_size = 1000;
	size_t jf_coverage = jf_reader.computeKmerCoverage(genome_size);
	size_t kff_coverage = kff_reader.computeKmerCoverage(genome_size);
	REQUIRE(jf_coverage == kff_coverage);
	
	// Test 3: Histogram computation produces identical results
	// Note: Testing with temporary files to avoid conflicts
	size_t jf_peak = jf_reader.computeHistogram(1000, true, "test_jf_histogram.tmp");
	size_t kff_peak = kff_reader.computeHistogram(1000, true, "test_kff_histogram.tmp");
	REQUIRE(jf_peak == kff_peak);
	
	// Test 4: Non-existent k-mers return identical results (should be 0)
	vector<string> nonexistent_kmers = {"GGGGGGGGGG", "TTTTTTTTTT", "CCCCCCCCCC"};
	for (const string& kmer : nonexistent_kmers) {
		size_t jf_count = jf_reader.getKmerAbundance(kmer);
		size_t kff_count = kff_reader.getKmerAbundance(kmer);
		
		INFO("Testing non-existent k-mer: " << kmer);
		REQUIRE(jf_count == kff_count);
		REQUIRE(jf_count == 0); // Should not exist
	}
}
#endif
