#pragma once
#ifndef __DEBUG_HELPER_H__
#define __DEBUG_HELPER_H__

#include <iostream>
#include <filesystem>
namespace debug_helper
{
	void print_bsr(const std::filesystem::path filename, const int& dim, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols);
}
#endif // __DEBUG_HELPER_H__
