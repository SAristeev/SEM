#include "dynamic_loads.h"
#include "num_types.h"
#include "solver.h"
#include "gll.h"
#include "parametric_hex.h"
#include "debug_helper.h"
#include "export2vtk.h"

#include <iostream>
#include <filesystem>
#include <cassert>

#include <mkl.h>
#include <mkl_spblas.h>

namespace solver{

	void buildFullGlobalMatrixStruct(const UnstructedMesh& mesh, std::vector<int>& rows, std::vector<int>& cols)
	{

		const std::vector<int>& shifts = mesh.elem_shifts;
		size_t elems_size = mesh.elem_type.size();
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
		gll::shape& shape_funcs = gll::shape::getInstance();
		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
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

		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++) {
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
						shape_funcs.dl(order, shape_index[0], point_index[0]) *
						shape_funcs.l(order, shape_index[1], point_index[1]) *
						shape_funcs.l(order, shape_index[2], point_index[2]);
					
					d_Nd_eta[jd] =
						shape_funcs.l(order, shape_index[0], point_index[0]) *
						shape_funcs.dl(order, shape_index[1], point_index[1]) *
						shape_funcs.l(order, shape_index[2], point_index[2]);

					d_Nd_zeta[jd] =
						shape_funcs.l(order, shape_index[0], point_index[0]) *
						shape_funcs.l(order, shape_index[1], point_index[1]) *
						shape_funcs.dl(order, shape_index[2], point_index[2]);
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

		
		/*Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
		auto rank = lu_decomp.rank();
		std::cout << "elem_id = " << elem_id << " rank = " << rank << " size = " << Bcols << std::endl;*/

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
		const UnstructedMesh& mesh = fcase.computational_mesh;
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
						const std::string& load_type = load.types[i].name;
						const std::vector<double>& params = load.types[i].param;
						double time_coeff = 1;
						if (load_type == "berlage")
						{
							time_coeff = berlage_load(0., params[0], params[1]);
						}
						F[dim * node + i] += load.data[i] * time_coeff;
					}
				}
			}
			if (load.name == "Pressure")
			{
				for (int p = 0; p < load.apply_to.size() / 2; p++) {
					int elem = 2 * p;
					int edge = 2 * p + 1;
					int shift1 = load.apply_to[edge];
					int shift2 = (load.apply_to[edge] + 1) % 4;
					int shift3 = (load.apply_to[edge] + 2) % 4;
					int shift4 = (load.apply_to[edge] + 3) % 4;

					int i = mesh.map_node_numeration.at(mesh.elems[mesh.elem_shifts[load.apply_to[elem] - 1] + shift1]);
					int j = mesh.map_node_numeration.at(mesh.elems[mesh.elem_shifts[load.apply_to[elem] - 1] + shift2]);
					int k = mesh.map_node_numeration.at(mesh.elems[mesh.elem_shifts[load.apply_to[elem] - 1] + shift3]);


					double nx = -(mesh.nodes[j][1] - mesh.nodes[i][1]);
					double ny = mesh.nodes[j][0] - mesh.nodes[i][0];

					if (nx * (mesh.nodes[k][0] - mesh.nodes[i][0]) + ny * (mesh.nodes[k][1] - mesh.nodes[i][1]) < 0)
					{
						double tmp = nx;
						nx = -ny;
						ny = tmp;
					}
					double len = std::sqrt(nx * nx + ny * ny);
					nx /= len;
					ny /= len;

					F[2 * i + 0] += load.data[0] * nx * len / 2;
					F[2 * i + 1] += load.data[0] * ny * len / 2;

					F[2 * j + 0] += load.data[0] * nx * len / 2;
					F[2 * j + 1] += load.data[0] * ny * len / 2;
				}
			}
		}
	}

	void updateLoads(const fc& fcase, std::vector<double>& F, double t)
	{
		const int& dim = fcase.dim;
		const UnstructedMesh& mesh = fcase.computational_mesh;
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
						const std::string& load_type = load.types[i].name;
						const std::vector<double>& params = load.types[i].param;
						double time_coeff = 1;
						if (load_type == "berlage") 
						{
							time_coeff = berlage_load(t, params[0], params[1]);
						}
						F[dim * node + i] += load.data[i] * time_coeff;
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
				int row = node;
				int n = fcase.computational_mesh.nodes.size();

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
	void LAE_solver(const int& blocksize, const std::vector<double>& A, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x) {

		int _iparm[64];
		void* _pt[64];
		{
			// Setup Pardiso control parameters
			for (int i = 0; i < 64; i++) {
				_iparm[i] = 0;
			}

			_iparm[0] = 1;  /* No solver default */
			_iparm[1] = 3; /* Fill-in reordering from METIS */ // !!! = 0
			/* Numbers of processors, value of OMP_NUM_THREADS */
			_iparm[2] = 1;
			_iparm[3] = 0; /* No iterative-direct algorithm */
			_iparm[4] = 0; /* No user fill-in reducing permutation */
			_iparm[5] = 0; /* If =0 then write solution only into x. If =1 then the RightHandSide-array will replaced by the solution*/
			_iparm[6] = 0; /* Not in use */
			_iparm[7] = 2; /* Max numbers of iterative refinement steps */
			_iparm[8] = 0; /* Not in use */
			_iparm[11] = 0; /* Not in use */
			if (0) { // not sym by default
				_iparm[9] = 8;
				_iparm[10] = 0; /* Disable scaling. Default for symmetric indefinite matrices. */
				_iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
			}
			else {
				_iparm[9] = 13;
				_iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
				_iparm[12] = 1; /*1!!! Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
			}
			//_iparm[12] = 1; /*1!!! Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
			_iparm[13] = 0; /* Output: Number of perturbed pivots */
			_iparm[14] = 0; /* Not in use */
			_iparm[15] = 0; /* Not in use */
			_iparm[16] = 0; /* Not in use */
			_iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
			_iparm[18] = -1; /* Output: Mflops for LU factorization */
			_iparm[19] = 0; /* Output: Numbers of CG Iterations */
			_iparm[26] = 1; /* Matrix Checker */
			if (1) // double by default
				_iparm[27] = 0;
			else
				_iparm[27] = 1;

			_iparm[59] = 0;

			_iparm[34] = 1; // zero-indexing
			if (blocksize == 1) {
				_iparm[36] = 0; // csr
			}
			else if (blocksize > 1) {
				_iparm[36] = blocksize; // bsr: block size
			}
			for (int i = 0; i < 64; i++) {
				_pt[i] = 0;
			}
		}
		MKL_INT n = rows.size() - 1;
		MKL_INT nnz = rows[n];

		const MKL_INT* h_RowsA = rows.data();
		const MKL_INT* h_ColsA = cols.data();
		const double* h_ValsA = A.data();

		const double* h_b = b.data();
		double* h_x = x.data();
		double ddum = 0.0;
		MKL_INT maxfct = 1;
		MKL_INT msglvl = 0;
		MKL_INT mnum = 1;
		MKL_INT mtype = 11;
		MKL_INT idum = 0;
		MKL_INT phase = 11;
		MKL_INT error = 0;

		//phase11

		pardiso(&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
			(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
			(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, &ddum, &ddum, (MKL_INT*)&error);


		//phase22
		phase = 22;
		pardiso((MKL_INT*)&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
			(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
			(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, &ddum, &ddum, (MKL_INT*)&error);

		//phase33
		phase = 33;
		pardiso((MKL_INT*)&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
			(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
			(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, (void*)h_b, (void*)h_x, (MKL_INT*)&error);

		//phase -1
		phase = -1;
		pardiso(&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
			(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
			(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, (void*)h_b, (void*)h_x, (MKL_INT*)&error);
		mkl_free_buffers();
	}

	void t_step(double dt, const int& blocksize, const std::vector<double>& M, const std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x, std::vector<double>& x_prev)
	{
		MKL_INT n = rows.size() - 1;
		MKL_INT nnz = rows[n];

		const MKL_INT* h_RowsK = rows.data();
		const MKL_INT* h_ColsK = cols.data();
		const double* h_ValsK = K.data();
		

		std::vector<double> tmp(x.size());
		const char trans = 'n';
		mkl_cspblas_dbsrgemv(&trans, &n, &blocksize, h_ValsK, h_RowsK, h_ColsK, x.data(), tmp.data());
		for (int i = 0; i < blocksize * n; i++)
		{
			tmp[i] = dt * dt * (b[i] - tmp[i]);
			x[i] = tmp[i] + x_prev[i];
		}
	}

	void solve(const fc& fcase, std::filesystem::path dir, std::string filename) 
	{
		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> K;
		buildFullGlobalMatrixStruct(fcase.computational_mesh, rows, cols);
		buildFullGlobalMatrix(fcase, K, rows, cols);
		std::vector<double> F;
		std::vector<double> x(fcase.dim * (rows.size() - 1));
		createLoads(fcase, F);
		applyconstraints(fcase, K, rows, cols, F);
		if (fcase.type == eDynamic) 
		{
			std::vector<double> x_prev = x;
			double t = 0;
			// TODO compute dt
			double dt = 5e-8;
			std::vector<double> time_steps;
			time_steps.push_back(0);
			// TODO compute M;
			std::vector<double> M(x.size());
			for (int i = 0; i < 1000; i++) 
			{
				std::cout << "step " << i << std::endl;
				updateLoads(fcase, F, t);
				t_step(dt, fcase.dim, M, K, rows, cols, 1, F, x, x_prev);
				post::export2vtk(fcase.computational_mesh, x, dir / std::string(filename + "_" + std::to_string(i) + ".vtu"));
				x_prev = x;
				t += dt;
				time_steps.push_back(t);			
			}
			post::collect_steps(dir, filename, time_steps);

		}
		else {
			LAE_solver(fcase.dim, K, rows, cols, 1, F, x);
			post::export2vtk(fcase.computational_mesh, x, dir / std::string(filename + ".fc"));
		}
	}
}
