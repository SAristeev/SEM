#include "fem.h"
#include "gll.h"
#include "num_types.h"
#include "dynamic_loads.h"
#include "gll.h"
#include "parametric_hex.h"
#include "debug_helper.h"
#include "export2vtk.h"

#include <mkl.h>
#include <cassert>

namespace solver 
{
	using namespace pre;
	
	

	static vtk_faces vtk_bound;

	void getMaterialMatrix(int dim, material_t mat, matd& D)
	{
		const int Ddim = 6;
		D = matd::Zero(Ddim, Ddim);

		// filling D
		const double& E = mat.E;
		const double& nu = mat.nu;
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
	}

	void computeJacobians(const fc& fcase, std::vector<std::vector<mat3>>& J, std::vector<std::vector<matd>>& B)
	{
		std::cout << "Compute Jacobians" << std::endl;


		//TODO: del clear
		J.clear();
		B.clear();

		gll::shape& shape_funcs = gll::shape::getInstance();
		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
		const material_t& material = fcase.materials[0];
		const int Ddim = 6;

		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++) {
			int order = mesh.order[elem_id];
			int nodes_per_edge = order - 1;
			const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];
			static const int Bcols = dim * nodes;

			B.emplace_back();
			B[elem_id].resize(nodes);

			J.emplace_back();
			J[elem_id].resize(nodes);

			std::vector<int> idx(nodes);
			for (int i = 0; i < nodes; i++)
			{
				idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
			}


			// J
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

				B[elem_id][id] = matd::Zero(Ddim, Bcols);
				matd Z = matd::Zero(Ddim, Bcols);
				J[elem_id][id] = mat3::Zero();

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

					J[elem_id][id](0, 0) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][0];
					J[elem_id][id](0, 1) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][1];
					J[elem_id][id](0, 2) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][2];

					J[elem_id][id](1, 0) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][0];
					J[elem_id][id](1, 1) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][1];
					J[elem_id][id](1, 2) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][2];

					J[elem_id][id](2, 0) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][0];
					J[elem_id][id](2, 1) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][1];
					J[elem_id][id](2, 2) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][2];
				}

				double detJ = J[elem_id][id].determinant();

				// TODO: debug logger (with levels)
				if (0) {
					std::cout << "mesh.elems[" << id << "] = " << idx[id] << " detJ = " << abs(detJ) << std::endl;
					if (abs(detJ) < 1e-14)
					{
						std::cout << J[elem_id][id] << std::endl;
					}
				}
				assert(abs(detJ) > 0);
				mat3 inv_J = J[elem_id][id].inverse();

				for (int jd = 0; jd < nodes; jd++)
				{
					d_Nd_x[jd] = inv_J(0, 0) * d_Nd_xi[jd] + inv_J(0, 1) * d_Nd_eta[jd] + inv_J(0, 2) * d_Nd_zeta[jd];
					d_Nd_y[jd] = inv_J(1, 0) * d_Nd_xi[jd] + inv_J(1, 1) * d_Nd_eta[jd] + inv_J(1, 2) * d_Nd_zeta[jd];
					d_Nd_z[jd] = inv_J(2, 0) * d_Nd_xi[jd] + inv_J(2, 1) * d_Nd_eta[jd] + inv_J(2, 2) * d_Nd_zeta[jd];
				}

				for (int jd = 0; jd < nodes; jd++)
				{
					B[elem_id][id](0, 3 * jd + 0) = d_Nd_x[jd];
					B[elem_id][id](1, 3 * jd + 0) = 0;
					B[elem_id][id](2, 3 * jd + 0) = 0;
					B[elem_id][id](3, 3 * jd + 0) = d_Nd_y[jd];
					B[elem_id][id](4, 3 * jd + 0) = 0;
					B[elem_id][id](5, 3 * jd + 0) = d_Nd_z[jd];

					B[elem_id][id](0, 3 * jd + 1) = 0;
					B[elem_id][id](1, 3 * jd + 1) = d_Nd_y[jd];
					B[elem_id][id](2, 3 * jd + 1) = 0;
					B[elem_id][id](3, 3 * jd + 1) = d_Nd_x[jd];
					B[elem_id][id](4, 3 * jd + 1) = d_Nd_z[jd];
					B[elem_id][id](5, 3 * jd + 1) = 0;

					B[elem_id][id](0, 3 * jd + 2) = 0;
					B[elem_id][id](1, 3 * jd + 2) = 0;
					B[elem_id][id](2, 3 * jd + 2) = d_Nd_z[jd];
					B[elem_id][id](3, 3 * jd + 2) = 0;
					B[elem_id][id](4, 3 * jd + 2) = d_Nd_y[jd];
					B[elem_id][id](5, 3 * jd + 2) = d_Nd_x[jd];

				}

			}
		}
	}

	void buildStiffnessMatrixStruct(const UnstructedMesh& mesh, std::vector<int>& rows, std::vector<int>& cols)
	{
		std::cout << "Builg Graph connectivity" << std::endl;
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

	void buildStiffnessMatrix(const fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<std::vector<mat3>>& J, const std::vector<std::vector<matd>>& B) {
		std::cout << "Building Stiffness matrix"<< std::endl;

		gll::shape& shape_funcs = gll::shape::getInstance();
		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
		const material_t& material = fcase.materials[0];

		size_t n = rows.size() - 1;
		size_t nnz = rows[n];

		K.resize(dim * dim * nnz);
		std::fill(K.begin(), K.end(), 0);

		const int Ddim = 6;
		matd D;
		getMaterialMatrix(dim, material, D);

		CBLAS_LAYOUT layout = CblasColMajor;
		CBLAS_TRANSPOSE nontrans = CblasNoTrans;
		CBLAS_TRANSPOSE trans = CblasTrans;
		const double beta = 1.0;

		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++) {
			int order = mesh.order[elem_id];
			int nodes_per_edge = order - 1;
			const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];
			static const int Bcols = dim * nodes;
			
			matd A = matd::Zero(Bcols, Bcols);
			matd Z = matd::Zero(Ddim, Bcols);

			std::vector<int> idx(nodes);
			for (int i = 0; i < nodes; i++)
			{
				idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
			}


			for (int id = 0; id < nodes; id++)
			{
				std::array<int, 3> point_index = pre::get_local_index(order, id);
				double weight = gll::weights[order - 1][point_index[0]] *
					gll::weights[order - 1][point_index[1]] *
					gll::weights[order - 1][point_index[2]];

				
				double detJ = J[elem_id][id].determinant();
				assert(abs(detJ) > 0);

				cblas_dgemm(layout, nontrans, nontrans, Ddim, Bcols, Ddim, std::abs(detJ), D.data(), Ddim, B[elem_id][id].data(), Ddim, 0, Z.data(), Ddim);
				cblas_dgemm(layout, trans, nontrans, Bcols, Bcols, Ddim, weight, B[elem_id][id].data(), Ddim, Z.data(), Ddim, 1, A.data(), Bcols);
			}

		
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

	void buildMassMatrix(const fc& fcase, std::vector<double>& M, const std::vector<std::vector<mat3>>& J, const std::vector<std::vector<matd>>& B) {
		std::cout << "Building Mass matrix" << std::endl;
		gll::shape& shape_funcs = gll::shape::getInstance();

		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
		const material_t& material = fcase.materials[0];

		M.resize(dim * mesh.nodes.size());
		std::fill(M.begin(), M.end(), 0);

		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++)
		{
			int order = mesh.order[elem_id];
			int nodes_per_edge = order - 1;
			const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];

			std::vector<int> idx(nodes);
			for (int i = 0; i < nodes; i++)
			{
				idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
			}

			for (int id = 0; id < nodes; id++)
			{
				std::array<int, 3> point_index = pre::get_local_index(order, id);
				
				double detJ = J[elem_id][id].determinant();
				double weight = gll::weights[order - 1][point_index[0]] *
								gll::weights[order - 1][point_index[1]] *
								gll::weights[order - 1][point_index[2]];
				M[dim * idx[id] + 0] += material.density * weight * std::abs(detJ);
				M[dim * idx[id] + 1] += material.density * weight * std::abs(detJ);
				M[dim * idx[id] + 2] += material.density * weight * std::abs(detJ);
			}
		}
	}

	void buildResultantMatrix(const fc& fcase, std::vector<double>& C, const std::vector<std::vector<mat3>>& J, const std::vector<std::vector<matd>>& B) {

		std::cout << "Building Resultant matrix" << std::endl;
		gll::shape& shape_funcs = gll::shape::getInstance();

		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
		
		C.resize(mesh.nodes.size());
		std::fill(C.begin(), C.end(), 0);

		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++)
		{
			int order = mesh.order[elem_id];
			int nodes_per_edge = order - 1;
			const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];

			std::vector<int> idx(nodes);
			for (int i = 0; i < nodes; i++)
			{
				idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
			}

			for (int id = 0; id < nodes; id++)
			{
				std::array<int, 3> point_index = pre::get_local_index(order, id);

				double detJ = J[elem_id][id].determinant();
				double weight = gll::weights[order - 1][point_index[0]] *
								gll::weights[order - 1][point_index[1]] *
								gll::weights[order - 1][point_index[2]];

				C[idx[id]] += weight * std::abs(detJ);
			}
		}
	}


	void resultants(const fc& fcase, std::vector<double>& C, std::vector<std::vector<double>>& eps, std::vector<std::vector<double>>& sigma, const std::vector<double>& u, const std::vector<std::vector<types::mat3>>& J, const std::vector<std::vector<types::matd>>& B) {


		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
		const material_t& material = fcase.materials[0];

		gll::shape& shape_funcs = gll::shape::getInstance();

		const int Ddim = 6;
		matd D;
		getMaterialMatrix(dim, material, D);

		CBLAS_LAYOUT layout = CblasColMajor;
		CBLAS_TRANSPOSE nontrans = CblasNoTrans;
		CBLAS_TRANSPOSE trans = CblasTrans;
		for (int kd = 0; kd < 6; kd++)
		{
			std::fill(eps[kd].begin(), eps[kd].end(), 0);
			std::fill(sigma[kd].begin(), sigma[kd].end(), 0);
		}
		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++) {
			int order = mesh.order[elem_id];
			int nodes_per_edge = order - 1;
			const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];

			std::vector<double> u_loc(dim * nodes);
			std::vector<double> sigma_loc(Ddim * nodes); // xx yy zz xy yz zx
			std::vector<double> eps_loc(Ddim * nodes); // xx yy zz xy yz zx

			std::vector<int> idx(nodes);
			for (int i = 0; i < nodes; i++)
			{
				idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
				u_loc[dim * i + 0] = u[dim * idx[i] + 0];
				u_loc[dim * i + 1] = u[dim * idx[i] + 1];
				u_loc[dim * i + 2] = u[dim * idx[i] + 2];
			}

			for (int id = 0; id < nodes; id++)
			{
				std::array<int, 3> point_index = pre::get_local_index(order, id);
				double weight = 
					gll::weights[order - 1][point_index[0]] *
					gll::weights[order - 1][point_index[1]] *
					gll::weights[order - 1][point_index[2]];

				cblas_dgemv(layout, nontrans, Ddim, dim * nodes, 1, B[elem_id][id].data(), Ddim, u_loc.data(), 1, 0.0, eps_loc.data(), 1);
				cblas_dgemv(layout, nontrans, Ddim, Ddim, 1, D.data(), Ddim, eps_loc.data(), 1, 0.0, sigma_loc.data(), 1);

				double detJ = J[elem_id][id].determinant();

				for (int kd = 0; kd < 6; kd++) 
				{
					eps[kd][idx[id]] += weight * eps_loc[kd] * std::abs(detJ);
					sigma[kd][idx[id]] += weight * sigma_loc[kd] * std::abs(detJ);
				}
			}
		}

		for (int kd = 0; kd < 6; kd++)
		{
			for (int node_id = 0; node_id < mesh.nodes.size(); node_id++) 
			{
				eps[kd][node_id] /= C[node_id];
				sigma[kd][node_id] /= C[node_id];
			}
		}
	}

	// BC

	void createLoads(const fc& fcase, std::vector<double>& F)
	{
	
		gll::shape& shape_funcs = gll::shape::getInstance();
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
				for (int pair_id = 0; pair_id < load.apply_to.size() / 2; pair_id++)
				{
					int elem = load.apply_to[2 * pair_id];
					int face = pre::faces_fidesys2vtk[load.apply_to[2 * pair_id + 1]];
					//int face = load.apply_to[2 * pair_id + 1];
					int order = mesh.order[elem];
					double pressure_data = load.data[0];

					std::vector<int> face_idxs;
					pre::get_face(order, face, face_idxs);

					int nodes = face_idxs.size();
					const std::array<int, 4>& local_ids = vtk_bound.faces[face];
					int elem_shift = mesh.elem_shifts[elem];
					
					vec3 p0 = mesh.nodes[mesh.elems[elem_shift + local_ids[0]]];
					vec3 p1 = mesh.nodes[mesh.elems[elem_shift + local_ids[1]]];
					vec3 p2 = mesh.nodes[mesh.elems[elem_shift + local_ids[2]]];

					vec3 v1 = p1 - p0;
					vec3 v2 = p2 - p0;

					int sign = face % 2 == 0 ? -1 : 1;
					vec3 outer_norm = (sign * v1.cross(v2)).normalized();
					//vec3 inner_norm = -outer_norm;
					vec3 pressure = outer_norm * pressure_data;
					
					int real_id_n = 0;

					for(int dim = 0; dim < 3; dim++){
						if (std::abs(outer_norm[dim]) > 1e-16) 
						{
							real_id_n = dim;
						}
					}
					int real_id_t1 = (real_id_n + 1) % 3;
					int real_id_t2 = (real_id_n + 2) % 3;

					vecd local_force = vecd::Zero(dim * nodes);

					std::vector<int> idx(nodes);
					for (int i = 0; i < nodes; i++)
					{
						idx[i] = mesh.elems[face_idxs[i] + mesh.elem_shifts[elem]];
					}

					mat2 J;

					// | dx/d\xi	dy/d\xi   |
					// | dx/d\eta	dy/d\eta  |	

					std::vector<double> d_Nd_xi(nodes);
					std::vector<double> d_Nd_eta(nodes);

					std::vector<double> d_Nd_x(nodes);
					std::vector<double> d_Nd_y(nodes);

					int loc_id_n = 0;

					if (face == 0 || face == 1) 
					{
						loc_id_n = 1;
					}
					else if (face == 2 || face == 3)
					{
						loc_id_n = 0;
					}
					else if (face == 4 || face == 5)
					{
						loc_id_n = 2;
					}
					else
					{
						assert(0);
					}

					int loc_id_t1 = (loc_id_n + 1) % 3;
					int loc_id_t2 = (loc_id_n + 2) % 3;

					for (int id = 0; id < nodes; id++) 
					{
						J = mat2::Zero();

						std::fill(d_Nd_xi.begin(), d_Nd_xi.end(), 0);
						std::fill(d_Nd_eta.begin(), d_Nd_eta.end(), 0);
						std::fill(d_Nd_x.begin(), d_Nd_x.end(), 0);
						std::fill(d_Nd_y.begin(), d_Nd_y.end(), 0);

						std::array<int, 3> point_index = pre::get_local_index(order, face_idxs[id]);
						for (int jd = 0; jd < nodes; jd++)
						{
							std::array<int, 3> shape_index = pre::get_local_index(order, face_idxs[jd]);

							d_Nd_xi[jd] =
								shape_funcs.dl(order, shape_index[loc_id_t1], point_index[loc_id_t1]) *
								shape_funcs.l(order, shape_index[loc_id_t2], point_index[loc_id_t2]);

							d_Nd_eta[jd] =
								shape_funcs.l(order, shape_index[loc_id_t1], point_index[loc_id_t1]) *
								shape_funcs.dl(order, shape_index[loc_id_t2], point_index[loc_id_t2]);
						}

						// compute J
						for (int jd = 0; jd < nodes; jd++)
						{
							J(0, 0) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][real_id_t1];
							J(0, 1) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][real_id_t2];

							J(1, 0) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][real_id_t1];
							J(1, 1) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][real_id_t2];
						}

						double detJ = J.determinant();
						
						//std::cout << detJ << std::endl;
						assert(abs(detJ) > 0);

						double weight = gll::weights[order - 1][point_index[loc_id_t1]] *
							gll::weights[order - 1][point_index[loc_id_t2]];

						F[dim * idx[id] + 0] += weight * pressure[0] * std::abs(detJ);
						F[dim * idx[id] + 1] += weight * pressure[1] * std::abs(detJ);
						F[dim * idx[id] + 2] += weight * pressure[2] * std::abs(detJ);
						
					}
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

	

	void applyconstraints(const fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, std::vector<double>& u)
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
					u[dim * row + 0] = constraint.data[0];
				}
				if (constraint.flag[1])
				{
					u[dim * row + 1] = constraint.data[1];
				}
				if (constraint.flag[2])
				{
					u[dim * row + 2] = constraint.data[1];
				}
			}
		}
	}
	
	void update_absorption(const fc& fcase, std::vector<double>& F,const std::vector<double>& velocity, const std::vector<double>& acceleration, const std::vector<int>& load_cut)
	{
		gll::shape& shape_funcs = gll::shape::getInstance();
		const material_t& mat = fcase.materials[0];

		double alpha = std::sqrt((mat.lambda + 2 * mat.mu) / mat.density);
		double beta = std::sqrt(mat.mu / mat.density);
		
		const UnstructedMesh& mesh = fcase.computational_mesh;
		int dim = fcase.dim;
		for (auto& load: fcase.loads) 
		{
			if (load.name == "absorption") 
			{
				for (int pair_id = 0; pair_id < load.apply_to.size() / 2; pair_id++)
				{
					int elem = load.apply_to[2 * pair_id];
					int face = pre::faces_fidesys2vtk[load.apply_to[2 * pair_id + 1]];
					//int face = load.apply_to[2 * pair_id + 1];
					int order = mesh.order[elem];


					std::vector<int> face_idxs;
					pre::get_face(order, face, face_idxs);

					int nodes = face_idxs.size();
					const std::array<int, 4>& local_ids = vtk_bound.faces[face];
					int elem_shift = mesh.elem_shifts[elem];

					vec3 p0 = mesh.nodes[mesh.elems[elem_shift + local_ids[0]]];
					vec3 p1 = mesh.nodes[mesh.elems[elem_shift + local_ids[1]]];
					vec3 p2 = mesh.nodes[mesh.elems[elem_shift + local_ids[2]]];

					vec3 v1 = p1 - p0;
					vec3 v2 = p2 - p0;

					int sign = face % 2 == 0 ? -1 : 1;
					vec3 outer_norm = (sign * v1.cross(v2)).normalized();
					vec3 inner_norm = -outer_norm;

					

					int real_id_n = 0;

					for (int dim = 0; dim < 3; dim++) {
						if (std::abs(outer_norm[dim]) > 1e-16)
						{
							real_id_n = dim;
						}
					}
					int real_id_t1 = (real_id_n + 1) % 3;
					int real_id_t2 = (real_id_n + 2) % 3;

					vecd local_force = vecd::Zero(dim * nodes);

					std::vector<int> idx(nodes);
					for (int i = 0; i < nodes; i++)
					{
						idx[i] = mesh.elems[face_idxs[i] + mesh.elem_shifts[elem]];
					}

					mat2 J;

					// | dx/d\xi	dy/d\xi   |
					// | dx/d\eta	dy/d\eta  |	

					std::vector<double> d_Nd_xi(nodes);
					std::vector<double> d_Nd_eta(nodes);

					std::vector<double> d_Nd_x(nodes);
					std::vector<double> d_Nd_y(nodes);

					int loc_id_n = 0;

					if (face == 0 || face == 1)
					{
						loc_id_n = 1;
					}
					else if (face == 2 || face == 3)
					{
						loc_id_n = 0;
					}
					else if (face == 4 || face == 5)
					{
						loc_id_n = 2;
					}
					else
					{
						assert(0);
					}

					int loc_id_t1 = (loc_id_n + 1) % 3;
					int loc_id_t2 = (loc_id_n + 2) % 3;

					for (int id = 0; id < nodes; id++)
					{
						J = mat2::Zero();

						std::fill(d_Nd_xi.begin(), d_Nd_xi.end(), 0);
						std::fill(d_Nd_eta.begin(), d_Nd_eta.end(), 0);
						std::fill(d_Nd_x.begin(), d_Nd_x.end(), 0);
						std::fill(d_Nd_y.begin(), d_Nd_y.end(), 0);

						std::array<int, 3> point_index = pre::get_local_index(order, face_idxs[id]);
						for (int jd = 0; jd < nodes; jd++)
						{
							std::array<int, 3> shape_index = pre::get_local_index(order, face_idxs[jd]);

							d_Nd_xi[jd] =
								shape_funcs.dl(order, shape_index[loc_id_t1], point_index[loc_id_t1]) *
								shape_funcs.l(order, shape_index[loc_id_t2], point_index[loc_id_t2]);

							d_Nd_eta[jd] =
								shape_funcs.l(order, shape_index[loc_id_t1], point_index[loc_id_t1]) *
								shape_funcs.dl(order, shape_index[loc_id_t2], point_index[loc_id_t2]);
						}

						// compute J
						for (int jd = 0; jd < nodes; jd++)
						{
							J(0, 0) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][real_id_t1];
							J(0, 1) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][real_id_t2];

							J(1, 0) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][real_id_t1];
							J(1, 1) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][real_id_t2];
						}

						double detJ = J.determinant();

						//std::cout << detJ << std::endl;
						assert(abs(detJ) > 0);

						double weight = gll::weights[order - 1][point_index[loc_id_t1]] * gll::weights[order - 1][point_index[loc_id_t2]];
						
						vec3 v = { velocity[dim * idx[id] + 0], velocity[dim * idx[id] + 1], velocity[dim * idx[id] + 2] };
						vec3 vn = v.dot(outer_norm) * outer_norm;
						vec3 vt = v - vn;
						vec3 result = - mat.density * (alpha * vn + beta * vt);
						
						
						F[dim * idx[id] + 0] += weight * result[0] * std::abs(detJ);
						F[dim * idx[id] + 1] += weight * result[1] * std::abs(detJ);
						F[dim * idx[id] + 2] += weight * result[2] * std::abs(detJ);

					}

				}

				
			}
		}
	}

	void updateconstraints(const fc& fcase, std::vector<double>& u)
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
				if (constraint.flag[0])
				{
					u[dim * node + 0] = constraint.data[0];
				}
				if (constraint.flag[1])
				{
					u[dim * node + 1] = constraint.data[1];
				}
				if (constraint.flag[2])
				{
					u[dim * node + 2] = constraint.data[2];
				}
			}
		}
	}

	

	
	// Second order Symmetric
	void explicit_step(const double dt, const int dim,const std::vector<double>& M, const std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<double>& b,  std::vector<double>& u_prev, std::vector<double>& u, std::vector<double>& u_next, std::vector<double>& Au)
	{

		int n = rows.size() - 1;
		int nnz = rows[n];
		const char trans = 'n';

		mkl_cspblas_dbsrgemv(&trans, &n, &dim, K.data(), rows.data(), cols.data(), u.data(), Au.data());
#pragma omp parallel for
		for (int i = 0; i < dim * n; i++) {

			Au[i] = dt * dt * (b[i] - Au[i]) / M[i];
			u_next[i] = Au[i] + 2 * u[i] - u_prev[i];
		}

	}

	void compute_velocity(const double dt, const int dim, std::vector<double>& v, const std::vector<double>& u_prev, const std::vector<double>& u, const std::vector<double>& u_next)
	{
		assert(v.size() != 0);
#pragma omp parallel for
		for (int i = 0; i < v.size(); i++) {
			v[i] = (u_next[i] - u_prev[i]) / (2 * dt);
		}
	}
	void compute_acceleration(const double dt, const int dim, std::vector<double>& a, const std::vector<double>& u_prev, const std::vector<double>& u, const std::vector<double>& u_next)
	{
		assert(a.size() != 0);
#pragma omp parallel for
		for (int i = 0; i < a.size(); i++) {
			a[i] = (u_next[i] - 2 * u[i] + u_prev[i]) / (dt * dt);
		}
	}









}
