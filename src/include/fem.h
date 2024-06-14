#pragma once
#ifndef __FEM_H__
#define __FEM_H__

#include "fc.h"

#include <vector>

namespace solver {
	
	void computeJacobians(const pre::fc& fcase, std::vector<std::vector<types::mat3>>& J, std::vector<std::vector<types::matd>>& B);

	void buildStiffnessMatrixStruct(const pre::UnstructedMesh& mesh, std::vector<int>& rows, std::vector<int>& cols);
	void buildStiffnessMatrix(const pre::fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<std::vector<types::mat3>>& J, const std::vector<std::vector<types::matd>>& B);

	void buildMassMatrix(const pre::fc& fcase, std::vector<double>& M, const std::vector<std::vector<types::mat3>>& J, const std::vector<std::vector<types::matd>>& B);
	void buildResultantMatrix(const pre::fc& fcase, std::vector<double>& C, const std::vector<std::vector<types::mat3>>& J, const std::vector<std::vector<types::matd>>& B);
	void resultants(const pre::fc& fcase, std::vector<double>& C, std::vector<std::vector<double>>& eps, std::vector<std::vector<double>>& sigma, const std::vector<double>& u, const std::vector<std::vector<types::mat3>>& J, const std::vector<std::vector<types::matd>>& B);

	void createLoads(const pre::fc& fcase, std::vector<double>& F);
	void updateLoads(const pre::fc& fcase, std::vector<double>& F, double t);
	void applyconstraints(const pre::fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, std::vector<double>& F);
	void updateconstraints(const pre::fc& fcase, std::vector<double>& u);
	void update_absorption(const pre::fc& fcase, std::vector<double>& F, const std::vector<double>& v, const std::vector<double>& a, const std::vector<int>& load_cut);
	

	void explicit_step(const double dt, const int dim,
		const std::vector<double>& M, const std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols,
		const std::vector<double>& b, std::vector<double>& u_prev, std::vector<double>& u, std::vector<double>& u_next,
		std::vector<double>& Au);
	void compute_velocity(const double dt, const int dim, std::vector<double>& v, const std::vector<double>& u_prev, const std::vector<double>& u, const std::vector<double>& u_next);
	void compute_acceleration(const double dt, const int dim, std::vector<double>& a, const std::vector<double>& u_prev, const std::vector<double>& u, const std::vector<double>& u_next);
}

#endif __FEM_H__