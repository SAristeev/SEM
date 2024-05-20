#include <num_types.h>
#include <solver.h>
#include <gll.h>
#include <parametric_hex.h>
#include <debug_helper.h>

#include <iostream>
#include <cassert>

#include <mkl.h>

namespace solver{

	void buildFullGlobalMatrixStruct(const UnstructedMesh& mesh, std::vector<int>& rows, std::vector<int>& cols)
	{

		const std::vector<int>& shifts = mesh.elem_shifts;
		size_t elems_size = mesh.elemids.size();
		size_t n = mesh.nodes.size();

		std::set<std::pair<int, int>> coostruct;

		for (int elem = 0; elem < elems_size; elem++)
		{
			for (int nodei = shifts[elem]; nodei < shifts[elem + 1]; nodei++)
			{
				for (int nodej = shifts[elem]; nodej < shifts[elem + 1]; nodej++)
				{
					coostruct.insert({ mesh.elems[nodei], mesh.elems[nodej] });
				}
			}
		}

		cols.resize(coostruct.size());
		rows.resize(n + 1);

		int currow = 0;
		int curnnz = rows[0] = 0;

		for (auto x : coostruct)
		{
			// zero based indexing
			if (x.first != currow)
			{
				rows[++currow] = curnnz;
			}
			cols[curnnz] = x.second;
			curnnz++;
		}
		rows[++currow] = curnnz;
	}

	void buildFullGlobalMatrix(const fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols) {

		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.mesh;
		const material_t& material = fcase.materials[0];

		size_t n = rows.size() - 1;
		size_t nnz = rows[n];

		K.resize(dim * dim * nnz);
		//std::fill(K.begin(), K.end(), 0);

		const int Ddim = 6;
		matd D = matd::Zero(Ddim, Ddim);
		
		// filling D
		const double& E = material.E;
		const double& nu = material.nu;
		const double D_const = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

		D(0, 0) = 1;
		D(1, 1) = 1;
		D(2, 2) = 1;

		D(0, 1) = nu / (1 - nu);
		D(0, 2) = nu / (1 - nu);

		D(1, 0) = nu / (1 - nu);
		D(1, 2) = nu / (1 - nu);

		D(2, 0) = nu / (1 - nu);
		D(2, 1) = nu / (1 - nu);

		D(3, 3) = (1 - 2 * nu) / (2 * (1 - nu));
		D(4, 4) = (1 - 2 * nu) / (2 * (1 - nu));
		D(5, 5) = (1 - 2 * nu) / (2 * (1 - nu));

		D *= D_const;

		CBLAS_LAYOUT layout = CblasColMajor;
		CBLAS_TRANSPOSE nontrans = CblasNoTrans;
		CBLAS_TRANSPOSE trans = CblasTrans;
		const double beta = 1.0;

		for (int elem_id = 0; elem_id < mesh.elemids.size(); elem_id++) {
			std::cout << "-------------------------------- " << std::endl;
			std::cout << "elem_id = " << elem_id << std::endl;
			int order = mesh.order[elem_id];
			int nodes_per_edge = order - 1;
			const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];
			static const int Bcols = dim * nodes;

			matd B = matd::Zero(Ddim, Bcols);
			matd A = matd::Zero(Bcols, Bcols);
			matd Z = matd::Zero(Ddim, Bcols);

			std::vector<int> idx(nodes);
			for (int i = 0; i < nodes; i++)
			{
				idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
			}

			Eigen::Matrix3d J;

			// | dx/d\xi	dy/d\xi     dz/d\xi   |
			// | dx/d\eta	dy/d\eta	dz/d\eta  |	
			// | dx/d\zeta	dy/d\zeta	dz/d\zeta |

			std::vector<double> d_Nd_xi(nodes);
			std::vector<double> d_Nd_eta(nodes);
			std::vector<double> d_Nd_zeta(nodes);

			std::vector<double> d_Nd_x(nodes);
			std::vector<double> d_Nd_y(nodes);
			std::vector<double> d_Nd_z(nodes);

			for (int id = 0; id < nodes; id++)
			{

				matd B = matd::Zero(Ddim, Bcols);
				matd Z = matd::Zero(Ddim, Bcols);
				J = mat3::Zero();
				
				std::fill(d_Nd_xi.begin(), d_Nd_xi.end(), 0);
				std::fill(d_Nd_eta.begin(), d_Nd_eta.end(), 0);
				std::fill(d_Nd_zeta.begin(), d_Nd_zeta.end(), 0);
				std::fill(d_Nd_x.begin(), d_Nd_x.end(), 0);
				std::fill(d_Nd_y.begin(), d_Nd_y.end(), 0);
				std::fill(d_Nd_z.begin(), d_Nd_z.end(), 0);
				
 				std::array<int, 3> point_index = pre::get_local_index(order, id);
				for (int jd = 0; jd < nodes; jd++)
				{
					std::array<int, 3> shape_index = pre::get_local_index(order, jd);

					d_Nd_xi[jd] =
						gll::func.dl(order, shape_index[0], point_index[0]) *
						gll::func.l(order, shape_index[1], point_index[1]) *
						gll::func.l(order, shape_index[2], point_index[2]);
					
					d_Nd_eta[jd] =
						gll::func.l(order, shape_index[0], point_index[0]) *
						gll::func.dl(order, shape_index[1], point_index[1]) *
						gll::func.l(order, shape_index[2], point_index[2]);

					d_Nd_zeta[jd] =
						gll::func.l(order, shape_index[0], point_index[0]) *
						gll::func.l(order, shape_index[1], point_index[1]) *
						gll::func.dl(order, shape_index[2], point_index[2]);
				}

				// compute J
				for (int jd = 0; jd < nodes; jd++)
				{

					J(0, 0) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][0];
					J(0, 1) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][1];
					J(0, 2) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][2];

					J(1, 0) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][0];
					J(1, 1) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][1];
					J(1, 2) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][2];

					J(2, 0) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][0];
					J(2, 1) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][1];
					J(2, 2) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][2];
				}

				double detJ = J.determinant();
				if(1){			
					std::cout << "mesh.elems[" << id << "] = " << idx[id]  << " detJ = " << abs(detJ) << std::endl;
					if (abs(detJ) < 1e-14) 
					{
						std::cout << J << std::endl;
					}
				}
				assert(abs(detJ) > 0);
				mat3 inv_J = J.inverse();

				for (int jd = 0; jd < nodes; jd++)
				{
					d_Nd_x[jd] = inv_J(0, 0) * d_Nd_xi[jd] + inv_J(0, 1) * d_Nd_eta[jd] + inv_J(0, 2) * d_Nd_zeta[jd];
					d_Nd_y[jd] = inv_J(1, 0) * d_Nd_xi[jd] + inv_J(1, 1) * d_Nd_eta[jd] + inv_J(1, 2) * d_Nd_zeta[jd];
					d_Nd_z[jd] = inv_J(2, 0) * d_Nd_xi[jd] + inv_J(2, 1) * d_Nd_eta[jd] + inv_J(2, 2) * d_Nd_zeta[jd];
				}

				for (int jd = 0; jd < nodes; jd++)
				{
					B(0, 3 * jd + 0) = d_Nd_x[jd];
					B(1, 3 * jd + 0) = 0;
					B(2, 3 * jd + 0) = 0;
					B(3, 3 * jd + 0) = d_Nd_y[jd];
					B(4, 3 * jd + 0) = 0;
					B(5, 3 * jd + 0) = d_Nd_z[jd];

					B(0, 3 * jd + 1) = 0;
					B(1, 3 * jd + 1) = d_Nd_y[jd];
					B(2, 3 * jd + 1) = 0;
					B(3, 3 * jd + 1) = d_Nd_x[jd];
					B(4, 3 * jd + 1) = d_Nd_z[jd];
					B(5, 3 * jd + 1) = 0;

					B(0, 3 * jd + 2) = 0;
					B(1, 3 * jd + 2) = 0;
					B(2, 3 * jd + 2) = d_Nd_z[jd];
					B(3, 3 * jd + 2) = 0;
					B(4, 3 * jd + 2) = d_Nd_y[jd];
					B(5, 3 * jd + 2) = d_Nd_x[jd];

				}
				cblas_dgemm(layout, nontrans, nontrans, Ddim, Bcols, Ddim, std::abs(detJ), D.data(), Ddim, B.data(), Ddim, 0, Z.data(), Ddim);
			
				cblas_dgemm(layout, trans, nontrans, Bcols, Bcols, Ddim, 1.0, B.data(), Ddim, Z.data(), Ddim, 1, A.data(), Bcols);
			}			

		
		Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
		auto rank = lu_decomp.rank();
		std::cout << "elem_id = " << elem_id << " rank = " << rank << " size = " << Bcols << std::endl;

			for (int id = 0; id < nodes; id++)
			{
				std::vector<int> glob_index(nodes);
				for (int q = rows[idx[id]]; q < rows[idx[id] + 1]; q++)
				{
					for (int kd = 0; kd < nodes; kd++)
					{
						if (cols[q] == idx[kd]) glob_index[kd] = q;
					}
				}
				for (int kd = 0; kd < nodes; kd++)
				{
					K[dim * dim * glob_index[kd] + 0] += A(dim * id + 0, dim * kd + 0);
					K[dim * dim * glob_index[kd] + 1] += A(dim * id + 0, dim * kd + 1);
					K[dim * dim * glob_index[kd] + 2] += A(dim * id + 0, dim * kd + 2);

					K[dim * dim * glob_index[kd] + 3] += A(dim * id + 1, dim * kd + 0);
					K[dim * dim * glob_index[kd] + 4] += A(dim * id + 1, dim * kd + 1);
					K[dim * dim * glob_index[kd] + 5] += A(dim * id + 1, dim * kd + 2);

					K[dim * dim * glob_index[kd] + 6] += A(dim * id + 2, dim * kd + 0);
					K[dim * dim * glob_index[kd] + 7] += A(dim * id + 2, dim * kd + 1);
					K[dim * dim * glob_index[kd] + 8] += A(dim * id + 2, dim * kd + 2);
				}
			}
		}		
	}

	void createLoads(const fc& fcase, std::vector<double>& F)
	{
		const int& dim = fcase.dim;
		const UnstructedMesh& mesh = fcase.mesh;
		const material_t& material = fcase.materials[0];

		F.resize(dim * mesh.nodes.size());
		std::fill(F.begin(), F.end(), 0);
		for (const BC& load : fcase.loads)
		{
			if (load.name == "Force")
			{
				for (int node : load.apply_to)
				{
					for (int i = 0; i < 3; i++)
					{
						F[2 * mesh.map_node_numeration.at(node) + i] += load.data[i];
					}
				}
			}
		}
	}

	// only for 1 order now
	void applyconstraints(const fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, std::vector<double>& F)
	{
		int dim = fcase.dim;
		assert(dim == 3);
		for (auto& constraint : fcase.restraints)
		{
			int coordinate = -1;
			for (int i = 0; i < 6; i++) {
				if (i >= 3 && constraint.flag[i] != false)
				{
					std::cout << "Found unsupported pivot constraint. Ignore " << std::endl;
				}
			}

			for (int node : constraint.apply_to)
			{
				int row = fcase.mesh.map_node_numeration.at(node);
				int n = fcase.mesh.nodes.size();

				// clear column
				for (int j = 0; j < n; j++)
				{
					for (int i = rows[j]; i < rows[j + 1]; i++)
					{
						if (constraint.flag[0] && cols[i] == row)
						{
							K[dim * dim * i + 0] = 0;
							K[dim * dim * i + 3] = 0;
							K[dim * dim * i + 6] = 0;
						}
						if (constraint.flag[1] && cols[i] == row)
						{
							K[dim * dim * i + 1] = 0;
							K[dim * dim * i + 4] = 0;
							K[dim * dim * i + 7] = 0;
						}
						if (constraint.flag[2] && cols[i] == row)
						{
							K[dim * dim * i + 2] = 0;
							K[dim * dim * i + 5] = 0;
							K[dim * dim * i + 8] = 0;
						}
					}
				}

				// clear row
				for (int i = rows[row]; i < rows[row + 1]; i++)
				{
					if (constraint.flag[0])
					{
						if (cols[i] == row)
						{
							K[dim * dim * i + 0] = 1;
						}
						else
						{
							K[dim * dim * i + 0] = 0;
						}
						K[dim * dim * i + 1] = 0;
						K[dim * dim * i + 2] = 0;
					}
					if (constraint.flag[1])
					{
						K[dim * dim * i + 3] = 0;
						if (cols[i] == row)
						{
							K[dim * dim * i + 4] = 1;
						}
						else
						{
							K[dim * dim * i + 4] = 0;
						}
						K[dim * dim * i + 5] = 0;

					}
					if (constraint.flag[2])
					{
						K[dim * dim * i + 6] = 0;
						K[dim * dim * i + 7] = 0;
						if (cols[i] == row)
						{
							K[dim * dim * i + 8] = 1;
						}
						else
						{
							K[dim * dim * i + 8] = 0;
						}
					}
				}
				if (constraint.flag[0])
				{
					F[dim * row + 0] = constraint.data[0];
				}
				if (constraint.flag[1])
				{
					F[dim * row + 1] = constraint.data[1];
				}
				if (constraint.flag[2])
				{
					F[dim * row + 2] = constraint.data[1];
				}
			}
		}
	}


	void solve(const fc& fcase) 
	{
		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> K;
		buildFullGlobalMatrixStruct(fcase.mesh, rows, cols);
		buildFullGlobalMatrix(fcase, K, rows, cols);

		std::vector<double> F;
		createLoads(fcase, F);
		applyconstraints(fcase, K, rows, cols, F);
		debug::print_bsr("C:/WD/Octave/Kc.txt", 3, K, rows, cols);
		int breakpoint = 0;
	}
}
