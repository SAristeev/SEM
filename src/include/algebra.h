#pragma once
#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include <vector>
namespace algebra 
{
	void bsr2csr(const int& blocksize, const std::vector<double>& bsr_A, const std::vector<int>& bsr_rows, const std::vector<int>& bsr_cols,
		std::vector<double>& csr_A, std::vector<int>& csr_rows, std::vector<int>& csr_cols);
	void solve_pardiso(const int& blocksize, const std::vector<double>& A, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x);
	void solve_amgcl(const int& blocksize, const std::vector<double>& A, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x);
}

#endif __ALGEBRA_H__