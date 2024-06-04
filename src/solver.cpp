#include "solver.h"
#include "fem.h"
#include "export2vtk.h"
#include "algebra.h"

#include <iostream>
#include <filesystem>
#include <cassert>


namespace solver{
	double compute_dt(const fc& fcase) 
	{
		const UnstructedMesh& mesh = fcase.computational_mesh;
		double min = (mesh.nodes[1] - mesh.nodes[0]).norm();
		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++)
		{
			for (int point_i = mesh.elem_shifts[elem_id]; point_i < mesh.elem_shifts[elem_id + 1]; point_i++)
			{
				for (int point_j = mesh.elem_shifts[elem_id]; point_j < mesh.elem_shifts[elem_id + 1]; point_j++)
				{
					if (point_i != point_j) 
					{
						min = std::min(min, (mesh.nodes[mesh.elems[point_i]] - mesh.nodes[mesh.elems[point_j]]).norm());
					}
				}
			}
		}
		const material_t& material = fcase.materials[0];
		double alpha = std::sqrt((material.lambda + 2 * material.mu) / material.density);
		double beta = std::sqrt(material.mu / material.density);

		std::cout << "min: " << min << std::endl;
		// V_max * dt / dx = q
		
		// v_R - speed of Rayleigh's wave (faster than S-wave and P-wave)
		// v_R ~= (beta * 0.87 + 0.12 * alpha) 
		// dt = q * |dx| / v_R
		
		double v_R = (beta * 0.87 + 0.12 * alpha);
		//double v_R = beta * (0.862 + 1.14 * material.nu) / (1 + material.nu);
		double v_max = std::max(std::max(alpha, beta), v_R);
		
		//return fcase.d_settings.courant * min * min / v_max;
		return 1e-3;
	}
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

		std::vector<double> u(problem_size, 0);
		if (fcase.type == eDynamic) 
		{
			std::cout << "Start Dynamic" << std::endl;
			std::vector<double> u_prev(problem_size, 0);
			std::vector<double> u_next(problem_size, 0);
			std::vector<double> v(problem_size, 0);
			std::vector<double> a(problem_size, 0);
			// specific for Spectral elements
			std::vector<double> M(problem_size, 0);
			buildMassMatrix(fcase, M);
			//double max_M = *std::max_element(M.begin(), M.end());
			//double min_M = *std::min_element(M.begin(), M.end());
			//std::cout << "max_u: " << max_M << " min_u: " << min_M << std::endl;
			// TODO compute dt from Courant Number
			std::vector<double> time_steps;
			std::vector<int> int_steps;
			double t = 0;
			double dt = 1e-3;// 7.34246377e-05;// compute_dt(fcase) / 20;
			
			
			// save init
			post::export2vtk(fcase.computational_mesh, u, v, F, dir / std::string(filename + "_0.vtu"));
			std::cout << "dt = " << dt << std::endl; 

			// if explicit
			if(1){		
				//updateLoads(fcase, F, dt);
				std::vector<double> buf1(problem_size, 0);
				for (int i = 1; i < 5000/*fcase.d_settings.max_iters*/; i++)
				{	
					// TODO: logger
					std::cout << "Step " << i << " Time: " << t << std::endl;

					t += dt;
					
					
					

					explicit_step(dt, fcase.dim, M, K, rows, cols, F, u_prev, u, u_next, buf1);
					std::fill(v.begin(), v.end(), 0);
					compute_velocity(dt, fcase.dim, v, u_prev, u, u_next);
					compute_acceleration(dt, fcase.dim, a, u_prev, u, u_next);

					updateLoads(fcase, F, t);
					update_absorption(fcase, F, v, u);
					updateconstraints(fcase, F);
					updateconstraints(fcase, v);
					updateconstraints(fcase, u_next);
					updateconstraints(fcase, u);
					updateconstraints(fcase, u_prev);

					double max_F = *std::max_element(F.begin(), F.end());
					double min_F = *std::min_element(F.begin(), F.end());

					double max_u = *std::max_element(u.begin(), u.end());
					double min_u = *std::min_element(u.begin(), u.end());

					double max_v = *std::max_element(u_prev.begin(), u_prev.end());
					double min_v = *std::min_element(u_prev.begin(), u_prev.end());
					if (1) 
					{
						std::cout << "max_u: " << max_u << " min_u: " << min_u << std::endl;
						std::cout << "max_v: " << max_v << " min_v: " << min_v << std::endl;
						std::cout << "max_F: " << max_F << " min_F: " << min_F << std::endl;
					}
					if(i % 5 == 0){
						int_steps.push_back(i);
						time_steps.push_back(t);
						// TODO: resultants
						post::export2vtk(fcase.computational_mesh, u, v, F, dir / std::string(filename + "_" + std::to_string(i) + ".vtu"));
					}
					u_prev = u;
					u = u_next;
					
					if (t > 100 * fcase.d_settings.max_time) 
					{
						std::cout << "Max time reached" << std::endl;
						break;
					}
				}
			}
			else 
			{
				//TODO: implicit (Newmark scheme)
			}
			post::collect_steps(dir, filename, int_steps, time_steps);
		}
		else {
			// constrants (Dirichlet BC)
			applyconstraints(fcase, K, rows, cols, F);

			algebra::solve_pardiso(fcase.dim, K, rows, cols, 1, F, u);
			post::export2vtk(fcase.computational_mesh, u, u, F, dir / std::string(filename + ".vtu"));
		}
	}
}


