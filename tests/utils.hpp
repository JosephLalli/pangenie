#include <vector>
#include <string>

bool doubles_equal(double a, double b);

bool compare_vectors (std::vector<double>& v1, std::vector<double>& v2);

void parse_vcf_lines(std::string filename, std::vector<std::vector<std::string>>& lines);

#ifdef KFF_SUPPORT
bool kff_files_equivalent(const std::string& jf_file, const std::string& kff_file, size_t k);
void create_equivalent_kff(const std::string& jf_file, const std::string& kff_file, size_t k);
bool verify_format_equivalence(const std::string& jf_path, const std::string& kff_path, size_t k);
std::vector<std::string> get_test_kmers();
#endif