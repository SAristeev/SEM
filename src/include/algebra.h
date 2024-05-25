#pragma once
#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include <vector>
namespace algebra 
{
	void solve_pardiso(const int& blocksize, const std::vector<double>& A, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x);
}

#endif __ALGEBRA_H__