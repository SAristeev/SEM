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

		std::vector<double> u(problem_size, 0);
		if (fcase.type == eDynamic) 
		{
			std::vector<double> u_prev(problem_size, 0);
			// specific for Spectral elements
			std::vector<double> M(problem_size, 0);
			buildMassMatrix(fcase, M);
			double max_M = *std::max_element(M.begin(), M.end());
			double min_M = *std::min_element(M.begin(), M.end());
			std::cout << "max_u: " << max_M << " min_u: " << min_M << std::endl;
			// TODO compute dt from Courant Number
			std::vector<double> time_steps;
			double t = 0;
			double dt = 1.95104799e-06;
			int max_iters = 500;
			
			// save init
			post::export2vtk(fcase.computational_mesh, u, dir / std::string(filename + "_0.vtu"));

			// if explicit
			if(1){		
				//updateLoads(fcase, F, dt);
				std::vector<double> buf1(problem_size, 0);
				std::vector<double> buf2(problem_size, 0);
				for (int i = 1; i < max_iters; i++)
				{	
					// TODO: logger
					std::cout << "Step " << i << " Time: " << t << std::endl;

					t += dt;
					
					
					updateLoads(fcase, F, t);

					explicit_step(dt, fcase.dim, M, K, rows, cols, F, u, u_prev, buf1, buf2);
					
					double max_F = *std::max_element(F.begin(), F.end());
					double min_F = *std::min_element(F.begin(), F.end());

					double max_u = *std::max_element(u.begin(), u.end());
					double min_u = *std::min_element(u.begin(), u.end());

					double max_v = *std::max_element(u_prev.begin(), u_prev.end());
					double min_v = *std::min_element(u_prev.begin(), u_prev.end());

					std::cout << "max_u: " << max_u << " min_u: " << min_u << std::endl;
					std::cout << "max_v: " << max_v << " min_v: " << min_v << std::endl;
					std::cout << "max_F: " << max_F << " min_F: " << min_F << std::endl;
					if(i % 1 == 0){
						time_steps.push_back(t);
						// TODO: resultants
						post::export2vtk(fcase.computational_mesh, u, dir / std::string(filename + "_" + std::to_string(i) + ".vtu"));
					}
				}
			}
			else 
			{
				//TODO: implicit (Newmark scheme)
			}
			post::collect_steps(dir, filename, time_steps);
		}
		else {
			algebra::solve_pardiso(fcase.dim, K, rows, cols, 1, F, u);
			post::export2vtk(fcase.computational_mesh, u, dir / std::string(filename + ".vtu"));
		}
	}
}


