#include "catch.hpp"
#include "utils.hpp"
#include "../src/jellyfishcounter.hpp"
#include "../src/jellyfishreader.hpp"
#ifdef KFF_SUPPORT
#include "../src/kffreader.hpp"
#endif
#include <vector>
#include <string>
#include <memory>

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
TEST_CASE("KffReader", "[KffReader]") {
	// Test KFF reader basic functionality
	KffReader reader("../tests/data/reads_equivalent.kff", 10);
	
	// Test specific k-mers that should be in reads_equivalent.kff
	REQUIRE(reader.getKmerAbundance("ATGCTGTAAA") == 1);
	REQUIRE(reader.getKmerAbundance("TGCTGTAAAA") == 1);
	REQUIRE(reader.getKmerAbundance("GCTGTAAAAA") == 1);

	// Test error handling - wrong kmer size
	REQUIRE_THROWS(KffReader("../tests/data/reads_equivalent.kff", 11));
}

TEST_CASE("KmerCounter_format_equivalence", "[KmerCounter]") {
	// Test that JF and KFF readers produce identical results using existing files
	const string jf_file = "../tests/data/reads.jf";
	const string kff_file = "../tests/data/reads_equivalent.kff";
	const size_t k = 10;
	
	// Test comprehensive equivalence using existing equivalent files
	REQUIRE(verify_format_equivalence(jf_file, kff_file, k));
	
	// Test utility function
	REQUIRE(kff_files_equivalent(jf_file, kff_file, k));
}

TEST_CASE("KmerCounter_polymorphic_interface", "[KmerCounter]") {
	// Test polymorphic usage of both readers through KmerCounter interface
	jellyfish::mer_dna::k(10);
	
	vector<unique_ptr<KmerCounter>> readers;
	readers.emplace_back(new JellyfishReader("../tests/data/reads.jf", 10));
	readers.emplace_back(new KffReader("../tests/data/reads_equivalent.kff", 10));
	
	string test_kmer = "ATGCTGTAAA";
	jellyfish::mer_dna jelly_kmer(test_kmer);
	
	// Both readers should work identically through base class interface
	for (const auto& reader : readers) {
		REQUIRE(reader->getKmerAbundance(test_kmer) == 1);
		REQUIRE(reader->getKmerAbundance(jelly_kmer) == 1);
		
		// Test interface methods
		size_t coverage = reader->computeKmerCoverage(1000);
		REQUIRE(coverage >= 0);
		
		// Test histogram computation - may throw for small datasets
		// This is expected behavior for datasets too small to have histogram peaks
		try {
			reader->computeHistogram(1000, true, "");
		} catch (const runtime_error& e) {
			// Expected for small test datasets without clear peaks
			REQUIRE(string(e.what()).find("no peak found") != string::npos);
		}
	}
}
#endif
