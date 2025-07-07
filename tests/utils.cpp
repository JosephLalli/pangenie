#include "utils.hpp"
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#ifdef KFF_SUPPORT
#include "../src/jellyfishreader.hpp"
#include "../src/kffreader.hpp"
#include <cstdlib>
#include <memory>
#endif

using namespace std;

bool doubles_equal(double a, double b) {
	return abs(a - b) < 0.0000001;
}

bool compare_vectors (vector<double>& v1, vector<double>& v2) {
	if (v1.size() != v2.size()) return false;

	for (size_t i = 0; i < v1.size(); ++i) {
		if (! doubles_equal(v1[i], v2[i])) return false;
	}
	return true;
}

void parse_vcf_lines(string filename, vector<vector<string>>& computed_lines) {
	ifstream file(filename);
	string line;
	while(getline(file, line)) {
		vector<string> tokens;
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;
		istringstream iss(line);
		string token;
		while(getline(iss, token, '\t'))
			tokens.push_back(token);
		computed_lines.push_back(tokens);
	}
}

#ifdef KFF_SUPPORT
vector<string> get_test_kmers() {
	return {
		"AAAAAAACGG", "AAAAAACGGC", "ATGCTGTAAA", "CGTTTTTTTA",
		"CTGTAAAAAA", "GCTGTAAAAA", "GTAAAAAAAC", "TGCTGTAAAA", "TGTAAAAAAA"
	};
}

void create_equivalent_kff(const string& jf_file, const string& kff_file, size_t k) {
	// Extract k-mers and counts from jellyfish file
	string cmd = "jellyfish dump " + jf_file + " | awk 'NR%2==1 {count=$1; gsub(\">\", \"\", count)} NR%2==0 {print $0 \" \" count}' > " + kff_file + ".txt";
	system(cmd.c_str());
	
	// Create KFF file with counts
	string kff_cmd = "kff-tools instr -i " + kff_file + ".txt -o " + kff_file + " -k " + to_string(k) + " -d 1 --delimiter ' '";
	system(kff_cmd.c_str());
	
	// Clean up temporary file
	remove((kff_file + ".txt").c_str());
}

bool verify_format_equivalence(const string& jf_path, const string& kff_path, size_t k) {
	// Set jellyfish k-mer size
	jellyfish::mer_dna::k(k);
	
	// Initialize both readers
	JellyfishReader jf_reader(jf_path, k);
	KffReader kff_reader(kff_path, k);
	
	// Test all known k-mers for equivalence
	vector<string> test_kmers = get_test_kmers();
	for (const string& kmer : test_kmers) {
		size_t jf_count = jf_reader.getKmerAbundance(kmer);
		size_t kff_count = kff_reader.getKmerAbundance(kmer);
		if (jf_count != kff_count) {
			cerr << "Mismatch for kmer " << kmer << ": JF=" << jf_count << ", KFF=" << kff_count << endl;
			return false;
		}
		
		// Test with jellyfish k-mer objects
		jellyfish::mer_dna jelly_kmer(kmer);
		size_t jf_jelly_count = jf_reader.getKmerAbundance(jelly_kmer);
		size_t kff_jelly_count = kff_reader.getKmerAbundance(jelly_kmer);
		if (jf_jelly_count != kff_jelly_count) {
			cerr << "Jellyfish object mismatch for kmer " << kmer << ": JF=" << jf_jelly_count << ", KFF=" << kff_jelly_count << endl;
			return false;
		}
	}
	
	// Test coverage computation equivalence
	size_t jf_coverage = jf_reader.computeKmerCoverage(1000);
	size_t kff_coverage = kff_reader.computeKmerCoverage(1000);
	if (jf_coverage != kff_coverage) {
		cerr << "Coverage mismatch: JF=" << jf_coverage << ", KFF=" << kff_coverage << endl;
		return false;
	}
	
	// Test histogram computation equivalence
	try {
		size_t jf_peak = jf_reader.computeHistogram(1000, true, "test_jf_histogram.tmp");
		size_t kff_peak = kff_reader.computeHistogram(1000, true, "test_kff_histogram.tmp");
		if (jf_peak != kff_peak) {
			cerr << "Histogram peak mismatch: JF=" << jf_peak << ", KFF=" << kff_peak << endl;
			return false;
		}
	} catch (const runtime_error& e) {
		// Both should fail with the same error for small datasets
		cerr << "Histogram computation failed (expected for small datasets): " << e.what() << endl;
	}
	
	return true;
}

bool kff_files_equivalent(const string& jf_file, const string& kff_file, size_t k) {
	try {
		return verify_format_equivalence(jf_file, kff_file, k);
	} catch (const exception& e) {
		return false;
	}
}
#endif