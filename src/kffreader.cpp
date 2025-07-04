#include "kffreader.hpp"

#ifdef KFF_SUPPORT

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <sstream>
#include "histogram.hpp"

using namespace std;

KffReader::KffReader(string kff_file, size_t kmersize) 
    : filename(kff_file), kmer_size(kmersize), cache_populated(false) {
    
    try {
        kff_reader = std::make_unique<Kff_reader>(kff_file);
        
        // Verify k-mer size matches
        if (kff_reader->get_var("k") != kmersize) {
            ostringstream oss;
            oss << "KffReader::KffReader: given kmer size (" << kmersize  
                << ") does not match KFF file kmer size (" << kff_reader->get_var("k") << ")." << endl;
            throw runtime_error(oss.str());
        }
        
        // Set jellyfish k-mer size for compatibility
        jellyfish::mer_dna::k(kmersize);
        
    } catch (const exception& e) {
        ostringstream oss;
        oss << "KffReader::KffReader: Failed to open KFF file '" << kff_file << "': " << e.what() << endl;
        throw runtime_error(oss.str());
    }
}

void KffReader::populate_cache() {
    if (cache_populated) {
        return;
    }
    
    kmer_cache.clear();
    
    // Reset reader to beginning
    kff_reader.reset();
    kff_reader = std::make_unique<Kff_reader>(filename);
    
    uint8_t* kmer_binary;
    uint8_t* data;
    
    // Iterate through all k-mers in the file
    while (kff_reader->has_next()) {
        kff_reader->next_kmer(kmer_binary, data);
        
        // Convert binary k-mer to string
        string kmer_str = kmer_to_string(kmer_binary, kmer_size);
        
        // Extract count from data (assuming first 8 bytes contain count)
        size_t count = 0;
        if (kff_reader->get_var("data_size") >= sizeof(size_t)) {
            count = *reinterpret_cast<size_t*>(data);
        } else if (kff_reader->get_var("data_size") >= sizeof(uint32_t)) {
            count = *reinterpret_cast<uint32_t*>(data);
        } else if (kff_reader->get_var("data_size") >= sizeof(uint16_t)) {
            count = *reinterpret_cast<uint16_t*>(data);
        } else if (kff_reader->get_var("data_size") >= sizeof(uint8_t)) {
            count = *reinterpret_cast<uint8_t*>(data);
        }
        
        kmer_cache[kmer_str] = count;
    }
    
    cache_populated = true;
}

string KffReader::kmer_to_string(uint8_t* kmer_binary, size_t k) {
    string result;
    result.reserve(k);
    
    // Convert binary k-mer to string using standard nucleotide encoding
    // Assuming 2-bit encoding: A=0, C=1, G=2, T=3
    for (size_t i = 0; i < k; i++) {
        uint8_t byte_idx = i / 4;
        uint8_t bit_offset = (i % 4) * 2;
        uint8_t nucleotide = (kmer_binary[byte_idx] >> bit_offset) & 0x3;
        
        switch (nucleotide) {
            case 0: result += 'A'; break;
            case 1: result += 'C'; break;
            case 2: result += 'G'; break;
            case 3: result += 'T'; break;
            default: result += 'N'; break;
        }
    }
    
    return result;
}

void KffReader::string_to_kmer(const string& kmer_str, uint8_t* kmer_binary) {
    size_t bytes_needed = (kmer_str.length() + 3) / 4;
    memset(kmer_binary, 0, bytes_needed);
    
    for (size_t i = 0; i < kmer_str.length(); i++) {
        uint8_t nucleotide = 0;
        switch (kmer_str[i]) {
            case 'A': case 'a': nucleotide = 0; break;
            case 'C': case 'c': nucleotide = 1; break;
            case 'G': case 'g': nucleotide = 2; break;
            case 'T': case 't': nucleotide = 3; break;
            default: nucleotide = 0; break; // Default to A for unknown nucleotides
        }
        
        uint8_t byte_idx = i / 4;
        uint8_t bit_offset = (i % 4) * 2;
        kmer_binary[byte_idx] |= (nucleotide << bit_offset);
    }
}

size_t KffReader::getKmerAbundance(string kmer) {
    if (!cache_populated) {
        populate_cache();
    }
    
    // Convert to uppercase for consistency
    transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
    
    // Try canonical form as well
    jellyfish::mer_dna jelly_kmer(kmer);
    jelly_kmer.canonicalize();
    string canonical_kmer = jelly_kmer.to_str();
    
    // Check both original and canonical forms
    auto it = kmer_cache.find(kmer);
    if (it != kmer_cache.end()) {
        return it->second;
    }
    
    it = kmer_cache.find(canonical_kmer);
    if (it != kmer_cache.end()) {
        return it->second;
    }
    
    return 0; // K-mer not found
}

size_t KffReader::getKmerAbundance(jellyfish::mer_dna jelly_kmer) {
    jelly_kmer.canonicalize();
    string kmer_str = jelly_kmer.to_str();
    return getKmerAbundance(kmer_str);
}

size_t KffReader::computeKmerCoverage(size_t genome_kmers) {
    if (!cache_populated) {
        populate_cache();
    }
    
    long double result = 0.0L;
    long double genome = 1.0L * genome_kmers;
    
    for (const auto& pair : kmer_cache) {
        long double count = 1.0L * pair.second;
        result += (count / genome);
    }
    
    return (size_t) ceil(result);
}

size_t KffReader::computeHistogram(size_t max_count, bool largest_peak, string filename) {
    if (!cache_populated) {
        populate_cache();
    }
    
    Histogram histogram(max_count);
    
    // Build histogram from cached k-mer counts
    for (const auto& pair : kmer_cache) {
        if (pair.second > 0) {
            histogram.add_value(pair.second);
        }
    }
    
    // Write histogram values to file
    if (filename != "") {
        histogram.write_to_file(filename);
    }
    
    // Smooth the histogram
    histogram.smooth_histogram();
    
    // Find peaks
    vector<size_t> peak_ids;
    vector<size_t> peak_values;
    histogram.find_peaks(peak_ids, peak_values);
    
    // Identify the largest and second largest (if it exists)
    if (peak_ids.size() == 0) {
        throw runtime_error("KffReader::computeHistogram: no peak found in kmer-count histogram.");
    }
    
    size_t kmer_coverage_estimate = -1;
    if (peak_ids.size() < 2) {
        cerr << "Histogram peak: " << peak_ids[0] << " (" << peak_values[0] << ")" << endl;
        kmer_coverage_estimate = peak_ids[0];
    } else {
        size_t largest, second, largest_id, second_id;
        if (peak_values[0] < peak_values[1]) {
            largest = peak_values[1];
            largest_id = peak_ids[1];
            second = peak_values[0];
            second_id = peak_ids[0];
        } else {
            largest = peak_values[0];
            largest_id = peak_ids[0];
            second = peak_values[1];
            second_id = peak_ids[1];
        }
        
        for (size_t i = 0; i < peak_values.size(); ++i) {
            if (peak_values[i] > largest) {
                second = largest;
                second_id = largest_id;
                largest = peak_values[i];
                largest_id = peak_ids[i];
            } else if ((peak_values[i] > second) && (peak_values[i] != largest)) {
                second = peak_values[i];
                second_id = peak_ids[i];
            }
        }
        
        cerr << "Histogram peaks: " << largest_id << " (" << largest << "), " << second_id << " (" << second << ")" << endl;
        if (largest_peak) {
            kmer_coverage_estimate = largest_id;
        } else {
            kmer_coverage_estimate = second_id;
        }
    }
    
    // Add expected abundance counts to end of hist file
    if (filename != "") {
        ofstream histofile;
        histofile.open(filename, ios::app);
        if (!histofile.good()) {
            stringstream ss;
            ss << "KffReader::computeHistogram: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
            throw runtime_error(ss.str());
        }
        histofile << "parameters\t" << kmer_coverage_estimate/2.0 << '\t' << kmer_coverage_estimate << endl;
        histofile.close();
    }
    
    return kmer_coverage_estimate;
}

KffReader::~KffReader() {
    kmer_cache.clear();
}

#endif // KFF_SUPPORT