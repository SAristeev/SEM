#include <debug_helper.h>
#include <fstream>

namespace debug {

	void print_bsr(const std::filesystem::path filename, const int& dim, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols) {
		std::ofstream file(filename);

		int n = rows.size() - 1;
		int nz = rows[n];

		file << n << " " << nz << " " << dim << std::endl;
		file << std::scientific << std::setprecision(20);
		for (size_t i = 0; i < n + 1; i++) { file << rows[i] << " "; } file << std::endl;
		for (size_t i = 0; i < nz; i++) { file << cols[i] << " "; } file << std::endl;
		for (size_t i = 0; i < dim * dim * nz; i++) { file << K[i] << " "; } file << std::endl;

		file.close();
	}
}