#ifndef KFFREADER_HPP
#define KFFREADER_HPP

#ifdef KFF_SUPPORT

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <jellyfish/mer_dna.hpp>
#include "kmercounter.hpp"
#include "kff_io.hpp"

/**
* Reads k-mer counts from KFF files using the KFF C++ API
**/

class KffReader : public KmerCounter {
public:
	/** 
	* @param kff_file path to the KFF file
	* @param kmersize expected k-mer size
	**/
	KffReader(std::string kff_file, size_t kmersize);
	
	/** get the abundance of given kmer (string) **/
	size_t getKmerAbundance(std::string kmer) override;

	/** get the abundance of given kmer (jellyfish kmer) **/
	size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) override;

	/** compute the kmer coverage relative to the number of kmers in the genome **/
	size_t computeKmerCoverage(size_t genome_kmers) override;

	/** computes kmer abundance histogram and returns the three highest peaks **/
	size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "") override;

	~KffReader();

private:
	/** name of input .kff file **/
	std::string filename;
	/** k-mer size **/
	size_t kmer_size;
	/** KFF file reader **/
	std::unique_ptr<Kff_reader> kff_reader;
	/** in-memory k-mer count cache for fast lookups **/
	std::map<std::string, size_t> kmer_cache;
	/** whether the cache has been populated **/
	bool cache_populated;
	
	/** populate the k-mer cache from KFF file **/
	void populate_cache();
	/** convert binary k-mer to string **/
	std::string kmer_to_string(uint8_t* kmer_binary, size_t k);
	/** convert string k-mer to binary **/
	void string_to_kmer(const std::string& kmer_str, uint8_t* kmer_binary);
};

#endif // KFF_SUPPORT
#endif // KFFREADER_HPP