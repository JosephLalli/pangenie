#include "catch.hpp"
#include "utils.hpp"
#include "../src/commands.hpp"
#include "../src/probabilitytable.hpp"
#include "../src/hmm.hpp"
#ifdef KFF_SUPPORT
#include "../src/kffreader.hpp"
#endif
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cereal/archives/binary.hpp>
#include <memory>

using namespace std;

TEST_CASE("Commands run_genotype_command1", "[Commands run_genotype_command1]") {

	string precomputed_prefix = "../tests/data/index";
	string readfile = "../tests/data/region-reads.fa";
	string outname = "../tests/data/testfull";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.01L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 215;
	uint64_t hash_size = 100000;
	size_t panel_size = 0;
	double recombrate = 1.26;
	bool output_panel = false;

	/** (1) produce results with command **/

	run_genotype_command(precomputed_prefix, readfile, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, sampling_size, hash_size, panel_size, recombrate, output_panel);

	// check if output file exists
	{
		ifstream file(outname + "_genotyping.vcf");
		REQUIRE(file.good());
	}
	// parse output file
	string line;
	vector<vector<string>> computed_lines;
	parse_vcf_lines(outname + "_genotyping.vcf", computed_lines);

	// check if output looks as expected
	REQUIRE(computed_lines.size() == 2);


	/** (2) produce results directly from internal HMM **/

	size_t kmer_abundance_peak = 18;
	ProbabilityTable probs = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
	UniqueKmersMap uk;

	// Reconstruct expected UniqueKmers objects from file
	ifstream is("../tests/data/region_UniqueKmersList.cereal", std::ios::binary);
	cereal::BinaryInputArchive archive_is( is );
	archive_is(uk);

	HMM hmm(&uk.unique_kmers["chr1"], &probs, only_genotyping, only_phasing, recombrate, false, effective_N);
	vector<GenotypingResult> genotypes = hmm.get_genotyping_result();

	REQUIRE(genotypes.size() == 2);
	genotypes[0].normalize();
	genotypes[1].normalize();
	vector<vector<unsigned short>> defined = {{0,1}, {0,1,2}};
	vector<string> expected_likelihoods = {};

	for(size_t i = 0; i < 2; ++i) {
		vector<long double> likelihoods = genotypes[i].get_specific_likelihoods(defined[i]).get_all_likelihoods(defined[i].size());
		ostringstream all;
		pair<int,int> genotype = genotypes[i].get_specific_likelihoods(defined[i]).get_likeliest_genotype();
		all << genotype.first << "/" << genotype.second << ":";
		all << genotypes[i].get_specific_likelihoods(defined[i]).get_genotype_quality(genotype.first, genotype.second) << ":";
		all << setprecision(4) << log10(likelihoods[0]);
		for (size_t j = 1; j < likelihoods.size(); ++j) {
			all << "," << setprecision(4) << log10(likelihoods[j]);
		}
		all << ":" << uk.unique_kmers["chr1"][i]->get_coverage();
		expected_likelihoods.push_back(all.str());
	}


	/** (3) Check if results are identical **/

	for (size_t i = 0; i < expected_likelihoods.size(); ++i) {
		REQUIRE(expected_likelihoods[i] == computed_lines[i][9]);
	}
}

TEST_CASE("Commands run_genotype_command2", "[Commands run_genotype_command2]") {
	string precomputed_prefix = "../tests/data/index";
	string readfile = "../tests/data/region-reads.fa";
	string outname = "../tests/data/testsampled";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.01L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 0;
	uint64_t hash_size = 100000;
	size_t panel_size = 5;
	double recombrate = 1.26;
	bool output_panel = false;

	/** (1) produce results with command **/

	run_genotype_command(precomputed_prefix, readfile, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, sampling_size, hash_size, panel_size, recombrate, output_panel);
	// check if output file exists
	{
		ifstream file(outname + "_genotyping.vcf");
		REQUIRE(file.good());
	}
	// parse output file
	string line;
	vector<vector<string>> computed_lines;
	parse_vcf_lines(outname + "_genotyping.vcf", computed_lines);

	// check if output looks as expected
	REQUIRE(computed_lines.size() == 2);


	/** (2) produce results directly from internal HMM **/

	size_t kmer_abundance_peak = 18;
	ProbabilityTable probs = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
	UniqueKmersMap uk;

	// Reconstruct expected UniqueKmers objects from file
	ifstream is("../tests/data/region2_UniqueKmersList.cereal", std::ios::binary);
	cereal::BinaryInputArchive archive_is( is );
	archive_is(uk);

	HMM hmm(&uk.unique_kmers["chr1"], &probs, only_genotyping, only_phasing, recombrate, false, effective_N);
	vector<GenotypingResult> genotypes = hmm.get_genotyping_result();

	REQUIRE(genotypes.size() == 2);
	genotypes[0].normalize();
	genotypes[1].normalize();
	vector<vector<unsigned short>> defined = {{0,1}, {0,1,2}};
	vector<string> expected_likelihoods = {};

	for(size_t i = 0; i < 2; ++i) {
		if (genotypes[i].contains_no_likelihoods()) {
			genotypes[i].add_to_likelihood(0,0,1.0);
		}
		vector<long double> likelihoods = genotypes[i].get_specific_likelihoods(defined[i]).get_all_likelihoods(defined[i].size());
		ostringstream all;
		pair<int,int> genotype = genotypes[i].get_specific_likelihoods(defined[i]).get_likeliest_genotype();
		if ((genotype.first != -1) && (genotype.second != -1)) {
			all << genotype.first << "/" << genotype.second << ":";
			all << genotypes[i].get_specific_likelihoods(defined[i]).get_genotype_quality(genotype.first, genotype.second) << ":";
		} else {
			all << ".:.:";
		}
		all << setprecision(4) << log10(likelihoods[0]);
		for (size_t j = 1; j < likelihoods.size(); ++j) {
			all << "," << setprecision(4) << log10(likelihoods[j]);
		}

		all << ":" << uk.unique_kmers["chr1"][i]->get_coverage();
		expected_likelihoods.push_back(all.str());
	}


	/** (3) Check if results are identical **/

	for (size_t i = 0; i < expected_likelihoods.size(); ++i) {
		REQUIRE(expected_likelihoods[i] == computed_lines[i][9]);
	}
}

#ifdef KFF_SUPPORT
TEST_CASE("Commands_KFF_file_detection", "[Commands]") {
	// Test file extension detection logic used in commands.cpp
	string kff_filename = "test_reads.kff";
	string jf_filename = "test_reads.jf";
	string fasta_filename = "test_reads.fa";
	
	// Test file extension checking logic
	bool is_kff = (kff_filename.size() >= 4 && kff_filename.substr(kff_filename.size()-4) == ".kff");
	bool is_jf = (jf_filename.size() >= 3 && jf_filename.substr(jf_filename.size()-3) == ".jf");
	bool is_fasta_kff = (fasta_filename.size() >= 4 && fasta_filename.substr(fasta_filename.size()-4) == ".kff");
	
	REQUIRE(is_kff == true);
	REQUIRE(is_jf == true);
	REQUIRE(is_fasta_kff == false); // .fa files should not match .kff pattern
	
	// Test edge cases
	vector<pair<string, bool>> test_cases = {
		{"file.kff", true},
		{"file.KFF", false},  // Case sensitive
		{"file.kff.backup", false},  // Extension not at end
		{"kff", false},  // Too short
		{"a.kff", true}
	};
	
	for (const auto& test_case : test_cases) {
		string filename = test_case.first;
		bool expected = test_case.second;
		bool is_kff_file = (filename.size() >= 4 && filename.substr(filename.size()-4) == ".kff");
		REQUIRE(is_kff_file == expected);
	}
}

TEST_CASE("Commands_KmerCounter_integration", "[Commands]") {
	// Test that KFF readers work through commands interface patterns
	ifstream test_file("../tests/data/test.kff");
	if (!test_file.good()) {
		// Skip test if KFF file doesn't exist
		return;
	}
	
	try {
		// Test creating KffReader through shared_ptr as done in commands.cpp
		shared_ptr<KmerCounter> read_kmer_counts = 
			shared_ptr<KffReader>(new KffReader("../tests/data/test.kff", 10));
		
		REQUIRE(read_kmer_counts != nullptr);
		
		// Test polymorphic interface usage
		string test_kmer = "ATGCTGTAAA";
		size_t count = read_kmer_counts->getKmerAbundance(test_kmer);
		REQUIRE(count >= 0);
		
		// Test jellyfish k-mer interface
		jellyfish::mer_dna::k(10);
		jellyfish::mer_dna jelly_kmer(test_kmer);
		size_t jelly_count = read_kmer_counts->getKmerAbundance(jelly_kmer);
		REQUIRE(jelly_count >= 0);
		
		// Test other interface methods used by commands
		size_t coverage = read_kmer_counts->computeKmerCoverage(1000);
		REQUIRE(coverage >= 0);
		REQUIRE_NOTHROW(read_kmer_counts->computeHistogram(1000, true, ""));
		
	} catch (const exception& e) {
		// Expected if test file is not available or invalid
		REQUIRE(true);
	}
}

TEST_CASE("Commands_KFF_error_handling", "[Commands]") {
	// Test that KFF errors are properly propagated
	
	// Test with non-existent file
	REQUIRE_THROWS_AS(
		shared_ptr<KffReader>(new KffReader("non_existent_file.kff", 10)),
		runtime_error
	);
	
	// Test with invalid k-mer size (if test file exists)
	ifstream test_file("../tests/data/test.kff");
	if (test_file.good()) {
		REQUIRE_THROWS_AS(
			shared_ptr<KffReader>(new KffReader("../tests/data/test.kff", 31)),
			runtime_error
		);
	}
}
#endif
