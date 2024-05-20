#pragma once
#ifndef __DEBUG_HELPER_H__
#define __DEBUG_HELPER_H__

#include <filesystem>
#include <tbb/tick_count.h>
namespace debug 
{
	// TODO: ADD MACRO MAGIC
	void print_bsr(const std::filesystem::path filename, const int& dim, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols);
}
#endif // __DEBUG_HELPER_H__
