#include "solver.h"
#include "fem.h"
#include "export2vtk.h"
#include "algebra.h"

#include <iostream>
#include <filesystem>
#include <cassert>


namespace solver{

	void start_problem(const fc& fcase, std::filesystem::path dir, std::string filename) 
	{
		int problem_size = fcase.dim * fcase.computational_mesh.nodes.size();
		// stiffness matrix
		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> K;
		buildStiffnessMatrixStruct(fcase.computational_mesh, rows, cols);
		buildStiffnessMatrix(fcase, K, rows, cols);
		
		// loads (Neumann BC)
		std::vector<double> F;
		createLoads(fcase, F);

		// constrants (Dirichlet BC)
		applyconstraints(fcase, K, rows, cols, F);

		std::vector<double> x(problem_size);
		if (fcase.type == eDynamic) 
		{
			std::vector<double> x_prev = x;
			// specific for Spectral elements
			std::vector<double> M(problem_size);
			// compute M

			// TODO compute dt from Courant Number
			std::vector<double> time_steps;
			double t = 0;
			double dt = 5e-8;
			int max_iters = 10;
			
			// save init
			post::export2vtk(fcase.computational_mesh, x, dir / std::string(filename + "_0.vtu"));

			// if explicit
			if(1){		
				
				std::vector<double> buf(problem_size);
				for (int i = 1; i < max_iters; i++)
				{
					// TODO: logger
					std::cout << "step " << i << std::endl;

					t += dt;
					time_steps.push_back(t);
					updateLoads(fcase, F, t);

					explicit_step(dt, fcase.dim, M, K, rows, cols, F, x, buf, x_prev);
					// TODO: resultants
					post::export2vtk(fcase.computational_mesh, x, dir / std::string(filename + "_" + std::to_string(i) + ".vtu"));
					x_prev = x;
				}
			}
			else 
			{
				//TODO:
				// implicit
				// Newmark scheme
			}
			post::collect_steps(dir, filename, time_steps);
		}
		else {
			algebra::solve_pardiso(fcase.dim, K, rows, cols, 1, F, x);
			post::export2vtk(fcase.computational_mesh, x, dir / std::string(filename + ".fc"));
		}
	}
}
