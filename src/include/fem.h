#pragma once
#ifndef __FEM_H__
#define __FEM_H__

#include "fc.h"

#include <vector>

namespace solver {
	void buildStiffnessMatrixStruct(const pre::UnstructedMesh& mesh, std::vector<int>& rows, std::vector<int>& cols);
	void buildStiffnessMatrix(const pre::fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols);
	void createLoads(const pre::fc& fcase, std::vector<double>& F);
	void updateLoads(const pre::fc& fcase, std::vector<double>& F, double t);
	void applyconstraints(const pre::fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, std::vector<double>& F);

	void explicit_step(const double dt, const int dim, const std::vector<double>& M, const std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<double>& b, std::vector<double>& x, std::vector<double>& buf, std::vector<double>& x_prev);
}

#endif __FEM_H__