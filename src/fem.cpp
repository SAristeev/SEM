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

	void buildStiffnessMatrix(const fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols) {
		std::cout << "Building Stiffnes matrix"<< std::endl;
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
				// TODO: debug logger (with levels)
				if (0) {
					std::cout << "mesh.elems[" << id << "] = " << idx[id] << " detJ = " << abs(detJ) << std::endl;
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
				double weight = gll::weights[order - 1][point_index[0]] *
								gll::weights[order - 1][point_index[1]] *
								gll::weights[order - 1][point_index[2]];
				cblas_dgemm(layout, trans, nontrans, Bcols, Bcols, Ddim, weight, B.data(), Ddim, Z.data(), Ddim, 1, A.data(), Bcols);
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

	void buildMassMatrix(const fc& fcase, std::vector<double>& M) {
		std::cout << "Building Mass matrix" << std::endl;
		gll::shape& shape_funcs = gll::shape::getInstance();
		const int& dim = fcase.dim;
		assert(dim == 3);
		const UnstructedMesh& mesh = fcase.computational_mesh;
		const material_t& material = fcase.materials[0];

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
				double weight = gll::weights[order - 1][point_index[0]] *
								gll::weights[order - 1][point_index[1]] *
								gll::weights[order - 1][point_index[2]];
				
				M[dim * idx[id] + 0] += material.density * weight * std::abs(detJ);
				M[dim * idx[id] + 1] += material.density * weight * std::abs(detJ);
				M[dim * idx[id] + 2] += material.density * weight * std::abs(detJ);
				
			}
		}
	}
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

					std::cout << "elem = " << elem;
					std::cout << " face = " << face;
					std::cout << " pressure = " << pressure_data;
					std::cout << std::endl;
					

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
					vec3 pressure = inner_norm * pressure_data;
					
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
					int face = load.apply_to[2 * pair_id + 1];
					int order = mesh.order[elem];
					std::vector<int> face_idxs;
					

					
					// vtk -> face_id
					// 0 - 1
					// 1 - 3
					// 2 - 0
					// 3 - 5
					// 4 - 2
					// 5 - 4

					// face_id -> vtk
					// 0 - 2
					// 1 - 0
					// 2 - 4
					// 3 - 1
					// 4 - 5
					// 5 - 3

					std::array<double, 6 > mask_fidesys = { 4, 1, 2, 3, 0, 5 };
					std::array<double, 6 > mask_2 =			{ 4, 1, 2, 1, 2, 5 };
					get_face(order, mask_2[face], face_idxs);
					int nodes = face_idxs.size();
					// 3 5 4 6 2 1 
					// 2 4 3 5 1 0
					// { 4, 1, 2, 3, 0, 5}; saved
					
					std::array<int, 4> face_4 = vtk_bound.faces[mask_fidesys[face]];

					int sign = 1;
					if (face == 4 || face == 3)
						sign = -1;

					vec3 p0 = mesh.nodes[mesh.elems[mesh.elem_shifts[elem] + face_4[0]]];
					vec3 p1 = mesh.nodes[mesh.elems[mesh.elem_shifts[elem] + face_4[1]]];
					vec3 p2 = mesh.nodes[mesh.elems[mesh.elem_shifts[elem] + face_4[2]]];

					vec3 t0 = mesh.nodes[mesh.elems[mesh.elem_shifts[elem] + face_idxs[0]]];
					vec3 t1 = mesh.nodes[mesh.elems[mesh.elem_shifts[elem] + face_idxs[1]]];
					vec3 t2 = mesh.nodes[mesh.elems[mesh.elem_shifts[elem] + face_idxs[2]]];

					vec3 v1 = p1 - p0;
					vec3 v2 = p2 - p0;

					vec3 n = sign * v1.cross(v2);
					n /= n.norm();

					// to integrate with gll
					std::vector<vec3> to_integrate(face_idxs.size());

					
					mat2 J;

					// | dx/d\xi	dy/d\xi  |
					// | dx/d\eta	dy/d\eta |
			

					
					std::vector<double> d_Nd_xi(nodes);
					std::vector<double> d_Nd_eta(nodes);

					std::vector<double> d_Nd_x(nodes);
					std::vector<double> d_Nd_y(nodes);
					
					std::array<bool, 3> mask_loc = { true, true, true };
					for (int id = 1; id < nodes; id++)
					{
						std::array<int, 3> ind = pre::get_local_index(order, face_idxs[id]);
						std::array<int, 3> jnd = pre::get_local_index(order, face_idxs[id]);

						mask_loc[0] = mask_loc[0] && ind[0] == jnd[0];
						mask_loc[1] = mask_loc[1] && ind[1] == jnd[1];
						mask_loc[2] = mask_loc[2] && ind[2] == jnd[2];
					}
					int loc_n = 0;
					int loc_t1 = 0;
					int loc_t2 = 0;
					for (int i = 0; i < 3; i++)
					{
						if (mask_loc[i])
						{
							loc_n = i;
							loc_t1 = (i + 1) % 3;
							loc_t2 = (i + 2) % 3;
						}
					}

					vec3 integral = vec3::Zero();
					for (int id = 0; id < nodes; id++)
					{
						int p_id = mesh.elems[mesh.elem_shifts[elem] + face_idxs[id]];
						std::array<int, 3> point_index = pre::get_local_index(order, face_idxs[id]);
						double size = (5 / 2);
						double detJ = size * size;// J.determinant();
						double weight = 1; // gll::weights[order - 1][point_index[loc_t1]] * gll::weights[order - 1][point_index[loc_t2]];

						vec3 v = { velocity[dim * p_id + 0], velocity[dim * p_id + 1], velocity[dim * p_id + 2] };
						vec3 a = { acceleration[dim * p_id + 0], acceleration[dim * p_id + 1], acceleration[dim * p_id + 2] };

						vec3 dvn_dt = v.dot(n) * v; // +v.dot(n) * a;
						vec3 dvr_dt = v - dvn_dt;

						//integral += mat.density * (alpha * dvn_dt + beta * dvr_dt) * detJ * weight;
						
						//vec3 dvn_dt = a.dot(n) * v + v.dot(n) * a;
						//vec3 dvr_dt = a - dvn_dt;
						integral += - mat.density *  (beta * dvn_dt + alpha * dvr_dt) * detJ * weight;
						//int doubled = 1;

						//for (const auto& constraint : fcase.restraints)
						//{
						//	const std::vector<int>& fixed_points = constraint.apply_to;
						//	//if (std::find(fixed_points.begin(), fixed_points.end(), p_id) != fixed_points.end()) 
						//	{
						//		//doubled *= 2;
						//		break;
						//	}
						//}
						//F[dim * p_id + 0] = doubled * to_integrate[id][0];
						//F[dim * p_id + 1] = doubled * to_integrate[id][1];
						//F[dim * p_id + 2] = doubled * to_integrate[id][2];
						
					}
					for (int id = 0; id < nodes; id++) 
					{
						int doubled = 1;
						int p_id = mesh.elems[mesh.elem_shifts[elem] + face_idxs[id]];
						for (const auto& constraint : fcase.restraints)
						{
							const std::vector<int>& fixed_points = constraint.apply_to;
							if (std::find(fixed_points.begin(), fixed_points.end(), p_id) != fixed_points.end()) 
							{
								doubled *= 2;
								break;
							}
							
						}
						
						if (p_id == 0 || p_id == 91) 
						{
							doubled = 4;
						}

						F[dim * p_id + 0] += doubled * integral[0];
						F[dim * p_id + 1] += doubled * integral[1];
						F[dim * p_id + 2] += doubled * integral[2];
					}
					int real_n = 0;
					int real_t1 = 0;
					int real_t2 = 0;
					for (int i = 0; i < 3; i++)
					{
						if (std::abs(n[i]) > 1e-10)
						{
							real_n = i;
							real_t1 = (i + 1) % 3;
							real_t2 = (i + 2) % 3;
						}
					}

					for (int id = 0; id < nodes; id++)
					{
						int p_id = mesh.elems[mesh.elem_shifts[elem] + face_idxs[id]];
						//F[dim * p_id + 0] += -acceleration[dim * p_id + 0];
						//F[dim * p_id + 1] += -acceleration[dim * p_id + 1];
						//F[dim * p_id + 2] += -acceleration[dim * p_id + 2];
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
	void explicit_step(const double dt, const int dim,
		const std::vector<double>& M, const std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols,
		const std::vector<double>& b, 
		std::vector<double>& u_prev, std::vector<double>& u, std::vector<double>& u_next,
		std::vector<double>& Au)
	{

		int n = rows.size() - 1;
		int nnz = rows[n];
		const char trans = 'n';

		mkl_cspblas_dbsrgemv(&trans, &n, &dim, K.data(), rows.data(), cols.data(), u.data(), Au.data());

		for (int i = 0; i < dim * n; i++) {
			if (std::abs(b[i]) > 1e-5) 
			{
				//std::cout << "b[i]" << b[i] << " Au[i] " << Au[i] << " i " << i << std::endl;
				//std::cout << "Au[i] / b[i] " << Au[i] / b[i] << " i " << i << std::endl;
			}
			Au[i] = dt * dt * (b[i] - Au[i]) / M[i];
			u_next[i] = Au[i] + 2 * u[i] - u_prev[i];
		}

	}
	void compute_velocity(const double dt, const int dim, std::vector<double>& v, const std::vector<double>& u_prev, const std::vector<double>& u, const std::vector<double>& u_next)
	{
		assert(v.size() != 0);
		for (int i = 0; i < v.size(); i++) {
			v[i] = (u_next[i] - u_prev[i]) / (2 * dt);
		}
	}
	void compute_acceleration(const double dt, const int dim, std::vector<double>& a, const std::vector<double>& u_prev, const std::vector<double>& u, const std::vector<double>& u_next)
	{
		assert(a.size() != 0);
		for (int i = 0; i < a.size(); i++) {
			a[i] = (u_next[i] - 2 * u[i] + u_prev[i]) / (dt * dt);
		}
	}


	//void resultants(const int& dim, material_t material, std::vector<double>& eps, std::vector<double>& sigma, const std::vector<double>& u, const UnstructedMesh& mesh, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols) {
	//	gll::shape& shape_funcs = gll::shape::getInstance();
	//	int blocksize = dim;
	//	size_t n = rows.size() - 1;
	//	size_t nnz = rows[n];
	//	size_t elems_size = mesh.elemids.size();

	//	std::vector<double> C(nnz, 0);
	//	std::vector<double> b(6 * n, 0);


	//	const int Ddim = 6;
	//	matd D = matd::Zero(Ddim, Ddim);

	//	// filling D
	//	const double& E = material.E;
	//	const double& nu = material.nu;
	//	const double D_const = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

	//	D(0, 0) = 1;
	//	D(1, 1) = 1;
	//	D(2, 2) = 1;

	//	D(0, 1) = nu / (1 - nu);
	//	D(0, 2) = nu / (1 - nu);

	//	D(1, 0) = nu / (1 - nu);
	//	D(1, 2) = nu / (1 - nu);

	//	D(2, 0) = nu / (1 - nu);
	//	D(2, 1) = nu / (1 - nu);

	//	D(3, 3) = (1 - 2 * nu) / (2 * (1 - nu));
	//	D(4, 4) = (1 - 2 * nu) / (2 * (1 - nu));
	//	D(5, 5) = (1 - 2 * nu) / (2 * (1 - nu));

	//	D *= D_const;


	//	CBLAS_LAYOUT layout = CblasColMajor;
	//	CBLAS_TRANSPOSE nontrans = CblasNoTrans;
	//	CBLAS_TRANSPOSE trans = CblasTrans;
	//	const double beta = 1.0;

	//	for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++) {
	//		int order = mesh.order[elem_id];
	//		int nodes_per_edge = order - 1;
	//		const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];
	//		static const int Bcols = dim * nodes;

	//		matd B = matd::Zero(Ddim, Bcols);
	//		matd A = matd::Zero(Bcols, Bcols);
	//		matd Z = matd::Zero(Ddim, Bcols);

	//		std::vector<int> idx(nodes);
	//		for (int i = 0; i < nodes; i++)
	//		{
	//			idx[i] = mesh.elems[i + mesh.elem_shifts[elem_id]];
	//		}

	//		Eigen::Matrix3d J;

	//		// | dx/d\xi	dy/d\xi     dz/d\xi   |
	//		// | dx/d\eta	dy/d\eta	dz/d\eta  |	
	//		// | dx/d\zeta	dy/d\zeta	dz/d\zeta |

	//		std::vector<double> d_Nd_xi(nodes);
	//		std::vector<double> d_Nd_eta(nodes);
	//		std::vector<double> d_Nd_zeta(nodes);

	//		std::vector<double> d_Nd_x(nodes);
	//		std::vector<double> d_Nd_y(nodes);
	//		std::vector<double> d_Nd_z(nodes);

	//		for (int id = 0; id < nodes; id++)
	//		{
	//			vecd u_loc = vecd::Zero(Bcols);
	//			vecd eps_loc = vecd::Zero(Bcols);
	//			vecd sigma_loc = vecd::Zero(Bcols);
	//			matd B = matd::Zero(Ddim, Bcols);
	//			matd Z = matd::Zero(Ddim, Bcols);
	//			J = mat3::Zero();

	//			std::fill(d_Nd_xi.begin(), d_Nd_xi.end(), 0);
	//			std::fill(d_Nd_eta.begin(), d_Nd_eta.end(), 0);
	//			std::fill(d_Nd_zeta.begin(), d_Nd_zeta.end(), 0);
	//			std::fill(d_Nd_x.begin(), d_Nd_x.end(), 0);
	//			std::fill(d_Nd_y.begin(), d_Nd_y.end(), 0);
	//			std::fill(d_Nd_z.begin(), d_Nd_z.end(), 0);

	//			std::array<int, 3> point_index = pre::get_local_index(order, id);
	//			for (int jd = 0; jd < nodes; jd++)
	//			{
	//				std::array<int, 3> shape_index = pre::get_local_index(order, jd);

	//				d_Nd_xi[jd] =
	//					shape_funcs.dl(order, shape_index[0], point_index[0]) *
	//					shape_funcs.l(order, shape_index[1], point_index[1]) *
	//					shape_funcs.l(order, shape_index[2], point_index[2]);

	//				d_Nd_eta[jd] =
	//					shape_funcs.l(order, shape_index[0], point_index[0]) *
	//					shape_funcs.dl(order, shape_index[1], point_index[1]) *
	//					shape_funcs.l(order, shape_index[2], point_index[2]);

	//				d_Nd_zeta[jd] =
	//					shape_funcs.l(order, shape_index[0], point_index[0]) *
	//					shape_funcs.l(order, shape_index[1], point_index[1]) *
	//					shape_funcs.dl(order, shape_index[2], point_index[2]);
	//			}

	//			// compute J
	//			for (int jd = 0; jd < nodes; jd++)
	//			{

	//				J(0, 0) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][0];
	//				J(0, 1) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][1];
	//				J(0, 2) += d_Nd_xi[jd] * mesh.nodes[idx[jd]][2];

	//				J(1, 0) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][0];
	//				J(1, 1) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][1];
	//				J(1, 2) += d_Nd_eta[jd] * mesh.nodes[idx[jd]][2];

	//				J(2, 0) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][0];
	//				J(2, 1) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][1];
	//				J(2, 2) += d_Nd_zeta[jd] * mesh.nodes[idx[jd]][2];
	//			}

	//			double detJ = J.determinant();
	//			assert(abs(detJ) > 0);
	//			mat3 inv_J = J.inverse();

	//			for (int jd = 0; jd < nodes; jd++)
	//			{
	//				d_Nd_x[jd] = inv_J(0, 0) * d_Nd_xi[jd] + inv_J(0, 1) * d_Nd_eta[jd] + inv_J(0, 2) * d_Nd_zeta[jd];
	//				d_Nd_y[jd] = inv_J(1, 0) * d_Nd_xi[jd] + inv_J(1, 1) * d_Nd_eta[jd] + inv_J(1, 2) * d_Nd_zeta[jd];
	//				d_Nd_z[jd] = inv_J(2, 0) * d_Nd_xi[jd] + inv_J(2, 1) * d_Nd_eta[jd] + inv_J(2, 2) * d_Nd_zeta[jd];
	//			}

	//			for (int jd = 0; jd < nodes; jd++)
	//			{
	//				B(0, 3 * jd + 0) = d_Nd_x[jd];
	//				B(1, 3 * jd + 0) = 0;
	//				B(2, 3 * jd + 0) = 0;
	//				B(3, 3 * jd + 0) = d_Nd_y[jd];
	//				B(4, 3 * jd + 0) = 0;
	//				B(5, 3 * jd + 0) = d_Nd_z[jd];

	//				B(0, 3 * jd + 1) = 0;
	//				B(1, 3 * jd + 1) = d_Nd_y[jd];
	//				B(2, 3 * jd + 1) = 0;
	//				B(3, 3 * jd + 1) = d_Nd_x[jd];
	//				B(4, 3 * jd + 1) = d_Nd_z[jd];
	//				B(5, 3 * jd + 1) = 0;

	//				B(0, 3 * jd + 2) = 0;
	//				B(1, 3 * jd + 2) = 0;
	//				B(2, 3 * jd + 2) = d_Nd_z[jd];
	//				B(3, 3 * jd + 2) = 0;
	//				B(4, 3 * jd + 2) = d_Nd_y[jd];
	//				B(5, 3 * jd + 2) = d_Nd_x[jd];

	//			}

	//			cblas_dgemv(layout, nontrans, Ddim, Bcols, std::abs(detJ), B.data(), Ddim, u.data(), 1, beta, eps_loc.data(), 1);
	//			cblas_dgemm(layout, nontrans, nontrans, Ddim, Bcols, Ddim, std::abs(detJ), D.data(), Ddim, B.data(), Ddim, 0, Z.data(), Ddim);
	//			eps_loc[2] /= 2;
	//			const double alpha = 1.0;
	//			const double beta = 0.0;
	//			cblas_dgemv(layout, nontrans, Ddim, Ddim, alpha, D.data(), Ddim, eps_loc.data(), 1, beta, sigma_loc.data(), 1);
	//			sigma_loc[2] *= 2;
	//			
	//		}



	//		
	//	}

	//	for (int e = 0; e < elems_size; e++) {
	//		std::fill(sigma_loc.begin(), sigma_loc.end(), 0.0);
	//		std::fill(eps_loc.begin(), eps_loc.end(), 0.0);

	//		int i = mesh.map_node_numeration.at(mesh.elems[mesh.nodes_per_elem[e] + 0]);
	//		int j = mesh.map_node_numeration.at(mesh.elems[mesh.nodes_per_elem[e] + 1]);
	//		int k = mesh.map_node_numeration.at(mesh.elems[mesh.nodes_per_elem[e] + 2]);


	//		double J = std::abs((mesh.nodes[j].x - mesh.nodes[i].x) * (mesh.nodes[k].y - mesh.nodes[i].y) -
	//			(mesh.nodes[k].x - mesh.nodes[i].x) * (mesh.nodes[j].y - mesh.nodes[i].y));

	//		
	//		u[0] = x[2 * i + 0];
	//		u[1] = x[2 * i + 1];
	//		u[2] = x[2 * j + 0];
	//		u[3] = x[2 * j + 1];
	//		u[4] = x[2 * k + 0];
	//		u[5] = x[2 * k + 1];

	//		CBLAS_LAYOUT layout = CblasColMajor;
	//		CBLAS_TRANSPOSE nontrans = CblasNoTrans;
	//		CBLAS_TRANSPOSE trans = CblasTrans;
	//		const double alpha = 1.0;
	//		const double beta = 0.0;
	//		cblas_dgemv(layout, nontrans, Brows, Bcols, 1 / (2 * S), B.data(), Brows, u.data(), 1, beta, eps_loc.data(), 1);
	//		eps_loc[2] /= 2;
	//		cblas_dgemv(layout, nontrans, Ddim, Ddim, alpha, D.data(), Ddim, eps_loc.data(), 1, beta, sigma_loc.data(), 1);
	//		sigma_loc[2] *= 2;
	//		int i1 = -1;
	//		int j1 = -1;
	//		int k1 = -1;
	//		for (int q = rows[i]; q < rows[i + 1]; q++) {
	//			//zero or one indexing
	//			if (cols[q] == i) {
	//				i1 = q;
	//			}
	//			if (cols[q] == j) {
	//				j1 = q;
	//			}
	//			if (cols[q] == k) {
	//				k1 = q;
	//			}
	//		}
	//		C[i1] += J / 12.0;
	//		C[j1] += J / 24.0;
	//		C[k1] += J / 24.0;

	//		b[i + 0 * n] += J * eps_loc[0] / 6.0;
	//		b[i + 1 * n] += J * eps_loc[1] / 6.0;
	//		b[i + 2 * n] += J * eps_loc[2] / 6.0;

	//		b[i + 3 * n] += J * sigma_loc[0] / 6.0;
	//		b[i + 4 * n] += J * sigma_loc[1] / 6.0;
	//		b[i + 5 * n] += J * sigma_loc[2] / 6.0;

	//		int i2 = -1;
	//		int j2 = -1;
	//		int k2 = -1;
	//		for (int q = rows[j]; q < rows[j + 1]; q++) {
	//			//zero or one indexing
	//			if (cols[q] == i) {
	//				i2 = q;
	//			}
	//			if (cols[q] == j) {
	//				j2 = q;
	//			}
	//			if (cols[q] == k) {
	//				k2 = q;
	//			}
	//		}
	//		C[i2] += J / 24.0;
	//		C[j2] += J / 12.0;
	//		C[k2] += J / 24.0;

	//		b[j + 0 * n] += J * eps_loc[0] / 6.0;
	//		b[j + 1 * n] += J * eps_loc[1] / 6.0;
	//		b[j + 2 * n] += J * eps_loc[2] / 6.0;

	//		b[j + 3 * n] += J * sigma_loc[0] / 6.0;
	//		b[j + 4 * n] += J * sigma_loc[1] / 6.0;
	//		b[j + 5 * n] += J * sigma_loc[2] / 6.0;

	//		int i3 = -1;
	//		int j3 = -1;
	//		int k3 = -1;
	//		for (int q = rows[k]; q < rows[k + 1]; q++) {
	//			//zero or one indexing
	//			if (cols[q] == i) {
	//				i3 = q;
	//			}
	//			if (cols[q] == j) {
	//				j3 = q;
	//			}
	//			if (cols[q] == k) {
	//				k3 = q;
	//			}
	//		}
	//		C[i3] += J / 24.0;
	//		C[j3] += J / 24.0;
	//		C[k3] += J / 12.0;

	//		b[k + 0 * n] += J * eps_loc[0] / 6.0;
	//		b[k + 1 * n] += J * eps_loc[1] / 6.0;
	//		b[k + 2 * n] += J * eps_loc[2] / 6.0;

	//		b[k + 3 * n] += J * sigma_loc[0] / 6.0;
	//		b[k + 4 * n] += J * sigma_loc[1] / 6.0;
	//		b[k + 5 * n] += J * sigma_loc[2] / 6.0;

	//	}



	//	const double* h_b = b.data();
	//	std::span<double> x_(reinterpret_cast<double*>(malloc(6 * n * sizeof(double))), 6 * n);
	//	double* h_x = x_.data();

	//	solve(1, C, rows, cols, 6, b, x_);

	//	eps = std::span(reinterpret_cast<double*>(malloc(3 * n * sizeof(double))), 3 * n);
	//	sigma = std::span(reinterpret_cast<double*>(malloc(3 * n * sizeof(double))), 3 * n);

	//	for (int i = 0; i < n; i++) {
	//		eps[3 * i + 0] = h_x[i + 0 * n];
	//		eps[3 * i + 1] = h_x[i + 1 * n];
	//		eps[3 * i + 2] = h_x[i + 2 * n];

	//		sigma[3 * i + 0] = h_x[i + 3 * n];
	//		sigma[3 * i + 1] = h_x[i + 4 * n];
	//		sigma[3 * i + 2] = h_x[i + 5 * n];
	//	}
	//	free(C.data());
	//	free(b.data());
	//	free(D.data());
	//	free(B.data());
	//	free(u.data());
	//	free(sigma_loc.data());

	//}

}
