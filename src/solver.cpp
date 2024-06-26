#include "debug_helper.h"
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

		
		// V_max * dt / dx = q
		
		// v_R - speed of Rayleigh's wave
		// v_R ~= (beta * 0.87 + 0.12 * alpha) 
		// dt = q * |dx| / v_R

		double v_R = (beta * 0.87 + 0.12 * alpha);
		double v_max = std::max(std::max(alpha, beta), v_R);
		
		return fcase.d_settings.courant * min / v_max;
		
	}
	void start_problem(const fc& fcase, std::filesystem::path dir, std::string filename) 
	{
		int dim = fcase.dim;
		int nodes_count = fcase.computational_mesh.nodes.size();
		

		std::vector<std::vector<matd>> B;
		std::vector<std::vector<mat3>> J;
		computeJacobians(fcase, J, B);

		// stiffness matrix
		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> K;
		
		buildStiffnessMatrixStruct(fcase.computational_mesh, rows, cols);
		buildStiffnessMatrix(fcase, K, rows, cols, J, B);
		
		std::vector<double> C;
		buildResultantMatrix(fcase, C, J, B);

		// loads (Neumann BC)
		std::vector<double> F;
		createLoads(fcase, F);
		
		std::vector<double> u(dim * nodes_count, 0);

		std::vector<std::vector<double>> eps(6); // xx yy zz xy yz zx
		std::vector<std::vector<double>> sigma(6); // xx yy zz xy yz zx

		for (int i = 0; i < 6; i++)
		{
			eps[i].resize(nodes_count);
			sigma[i].resize(nodes_count);
		}

		std::filesystem::create_directory(dir / "results");
		
		
		if (fcase.type == eDynamic) 
		{
			std::cout << "Start Dynamic" << std::endl;
			std::vector<double> u_prev(dim * nodes_count, 0);
			std::vector<double> u_next(dim * nodes_count, 0);
			std::vector<double> v(dim * nodes_count, 0);
			std::vector<double> a(dim * nodes_count, 0);

			
			// specific for Spectral elements
			std::vector<double> M;
			buildMassMatrix(fcase, M, J, B);

			// TODO compute dt from Courant Number
			std::vector<double> time_steps;
			std::vector<int> int_steps;
			double t = 0;
			//double dt = 5e-6;// compute_dt(fcase);
			double dt = 2.5e-6;
			double maxTime = 0.05;
			int max_iters = static_cast<int>(maxTime / dt);
			// save init
			
			post::export2vtk(fcase.computational_mesh, u, v, a, eps, sigma, F, dir / "results" / std::string(filename + "_0.vtu"));
			std::cout << "dt = " << dt << std::endl; 

			
			if(1){		
			
				std::vector<double> buf1(dim * nodes_count, 0);
				for (int i = 1; i <= max_iters + 1; i++)
				{	
					// TODO: logger
					std::cout << "Step " << i << " Time: " << t << std::endl;

					t += dt;
					
					updateLoads(fcase, F, t);
					//update_absorption(fcase, F, u, a, load_cut);

					explicit_step(dt, fcase.dim, M, K, rows, cols, F, u_prev, u, u_next, buf1);
					
					updateconstraints(fcase, u_next);
					updateconstraints(fcase, u);
					updateconstraints(fcase, u_prev);

					compute_velocity(dt, fcase.dim, v, u_prev, u, u_next);
					compute_acceleration(dt, fcase.dim, a, u_prev, u, u_next);
					
					

					double max_F = *std::max_element(F.begin(), F.end());
					double min_F = *std::min_element(F.begin(), F.end());

					double max_u = *std::max_element(u.begin(), u.end());
					double min_u = *std::min_element(u.begin(), u.end());

					double max_v = *std::max_element(u_prev.begin(), u_prev.end());
					double min_v = *std::min_element(u_prev.begin(), u_prev.end());
					
					double max_a = *std::max_element(a.begin(), a.end());
					double min_a = *std::min_element(a.begin(), a.end());
					if (1) 
					{
						std::cout << "max_u: " << max_u << " min_u: " << min_u << std::endl;
						std::cout << "max_v: " << max_v << " min_v: " << min_v << std::endl;
						std::cout << "max_a: " << max_a << " min_a: " << min_a << std::endl;
						std::cout << "max_F: " << max_F << " min_F: " << min_F << std::endl;
					}

					if (i % 100 == 0) {
						int_steps.push_back(i);
						time_steps.push_back(t);
						
						resultants(fcase, C, eps, sigma, u, J, B);
						post::export2vtk(fcase.computational_mesh, u, v, a, eps, sigma, F, dir / "results" / std::string(filename + "_" + std::to_string(i) + ".vtu"));
					}
					u_prev = u;
					u = u_next;
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
			applyconstraints(fcase, K, rows, cols, u);
			// stiffness matrix
			std::vector<int> csr_rows;
			std::vector<int> csr_cols;
			std::vector<double> csr_K;
			algebra::bsr2csr(fcase.dim, K, rows, cols, csr_K, csr_rows, csr_cols);
			//debug_helper::print_bsr("C:/WD/octave/bsr.txt", 3, K, rows, cols);
			//debug_helper::print_bsr("C:/WD/octave/csr.txt", 1, csr_K, csr_rows, csr_cols);
			algebra::solve_amgcl(1, csr_K, csr_rows, csr_cols, 1, F, u);

			resultants(fcase, C, eps, sigma, u, J, B);
			post::export2vtk(fcase.computational_mesh, u, u, u, eps, sigma, F, dir / "results" / std::string(filename + ".vtu"));
		}
	}
}


