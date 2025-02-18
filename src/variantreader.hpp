#ifndef VARIANT_READER_HPP
#define VARIANT_READER_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>
#include "fastareader.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"
#include "uniquekmers.hpp"
#include "dnasequence.hpp"
#include <functional>
#include <filesystem>

#include "cereal/access.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/map.hpp"
#include <cereal/archives/binary.hpp>


//std::vector<unsigned char> construct_index(std::vector<DnaSequence>& alleles, bool reference_added);
//std::vector<unsigned char> construct_index(std::vector<std::string>& alleles, bool reference_added);

template<class T>
std::vector<unsigned char> construct_index(std::vector<T>& alleles, bool reference_added) {
	size_t length = alleles.size();
	unsigned char offset = 0;
	if (reference_added) {
		assert(length > 0);
		length -= 1;
		offset += 1;
	}
	std::vector<unsigned char> index(length);
	std::iota(index.begin(), index.end(), 0);
	std::sort(index.begin(), index.end(), [&](unsigned char a, unsigned char b) { return alleles[a+offset] < alleles[b+offset]; });
	return index;
}

std::string hash_filenames(std::string reference, std::string vcf);

class VariantReader {
public:
    VariantReader() = default;
	VariantReader(std::string filename);
	VariantReader(std::string filename, std::string reference_filename, size_t kmer_size, bool add_reference, std::string sample = "sample");
	/**  writes all path segments (allele sequences + reference sequences in between)
	*    to the given file.
	**/
    FastaReader fasta_reader;
	size_t get_kmer_size() const;
	void write_path_segments(std::string filename) const;
	void get_chromosomes(std::vector<std::string>* result) const;
	size_t size_of(std::string chromosome) const;
	const Variant& get_variant(std::string chromosome, size_t index) const;
	const std::vector<Variant>& get_variants_on_chromosome(std::string chromosome) const;
	void write_outfiles(const Results &results, UniqueKmersMap &unique_kmers_list, std::string gt_outfile, std::string phase_outfile, bool only_genotyping, bool only_phasing, bool ignore_imputed);
	void write_genotypes_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, std::vector<UniqueKmers*>* unique_kmers, std::ofstream &gt_stream, bool ignore_imputed = false);
	void write_phasing_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, std::vector<UniqueKmers*>* unique_kmers, std::ofstream &phasing_stream, bool ignore_imputed = false);
	// size_t nr_of_genomic_kmers() const;
	size_t nr_of_paths() const;
	void get_left_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
	void get_right_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
    void Store(std::string filename) const;
    void Load(std::string name);
	void swap(VariantReader & var) noexcept;
	VariantReader & operator=(VariantReader var);
	std::string sample;

private:
    friend cereal::access;
	template<class Archive>
	void serialize(Archive & archive) {
    	archive(kmer_size, nr_paths, nr_variants, add_reference, sample, variants_per_chromosome, variant_ids); 
    }
	void open_genotyping_outfile(std::string outfile_name, std::ofstream &genotyping_outfile);
	void open_phasing_outfile(std::string outfile_name, std::ofstream &phasing_outfile);
	void close_genotyping_outfile(std::ofstream &gt_stream);
	void close_phasing_outfile(std::ofstream &phasing_outfile);

    std::string REF_VCF_HASH_NAME;
	size_t kmer_size;
	size_t nr_paths;
	size_t nr_variants;
	bool add_reference;
	bool genotyping_outfile_open;
	bool phasing_outfile_open;
	std::map< std::string, std::vector<Variant> > variants_per_chromosome;
	std::map< std::string, std::vector<std::vector<std::string>>> variant_ids;
	void add_variant_cluster(std::string& chromosome, std::vector<Variant>* cluster);
	void insert_ids(std::string& chromosome, std::vector<DnaSequence>& alleles, std::vector<std::string>& variant_ids, bool reference_added);
	std::string get_ids(std::string chromosome, std::vector<std::string>& alleles, size_t variant_index, bool reference_added);
};

#endif // VARIANT_READER_HPP
