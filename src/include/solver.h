#pragma once
#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "fc.h"
namespace solver {
	using namespace pre;
	void start_problem(const fc& fcase, std::filesystem::path dir, std::string filename);
}

#endif __SOLVER_H__