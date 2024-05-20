#include<solver.h>

#include <cassert>

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
			if (x.first == currow + 1)
			{
				rows[++currow] = curnnz;
			}
			cols[curnnz] = x.second;
			curnnz++;
		}
		rows[++currow] = curnnz;
	}

	//void buildFullGlobalMatrix(const fc& fcase, std::vector<double>& K, const std::vector<int>& rows, const std::vector<int>& cols) {

	//	const int& dim = fcase.dim;
	//	assert(dim == 2);
	//	const UnstructedMesh& mesh = fcase.mesh;
	//	const material_t& material = fcase.materials[0];

	//	int blocksize = dim;
	//	size_t n = rows.size() - 1;
	//	size_t nnz = rows[n];

	//	K.resize(blocksize * blocksize * nnz);
	//	std::fill(K.begin(), K.end(), 0);

	//	const int Ddim = 6;
	//	Eigen::Matrix<double, Ddim, Ddim> D = Eigen::Matrix<double, Ddim, Ddim>::Zero();

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

	//	for (int elem_id = 0; elem_id < mesh.elemids.size(); elem_id++) {
	//		int order = mesh.order[elem_id];
	//		int nodes_per_edge = order - 1;
	//		const int nodes = mesh.elem_shifts[elem_id + 1] - mesh.elem_shifts[elem_id];
	//		static const int Bcols = dim * nodes;

	//		Eigen::MatrixXd B = Eigen::MatrixXd::Zero(Ddim, Bcols);
	//		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Bcols, Bcols);
	//		Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(Ddim, Bcols);

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

	//			Eigen::MatrixXd B = Eigen::MatrixXd::Zero(Ddim, Bcols);
	//			Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(Ddim, Bcols);

	//			std::fill(d_Nd_xi.begin(), d_Nd_xi.end(), 0);
	//			std::fill(d_Nd_eta.begin(), d_Nd_eta.end(), 0);
	//			std::fill(d_Nd_zeta.begin(), d_Nd_zeta.end(), 0);
	//			std::fill(d_Nd_x.begin(), d_Nd_x.end(), 0);
	//			std::fill(d_Nd_y.begin(), d_Nd_y.end(), 0);
	//			std::fill(d_Nd_z.begin(), d_Nd_z.end(), 0);

	//			J = Eigen::Matrix3d::Zero();
	//			//std::fill(J.begin(), J.end(), 0);

	//			std::tuple<int, int, int> m_id = mesh.hex_loc_indxs.at(id);
	//			int i_xi = get<0>(m_id);
	//			int i_eta = get<1>(m_id);
	//			int i_zeta = get<2>(m_id);
	//			for (int jd = 0; jd < nodes; jd++)
	//			{
	//				std::tuple<int, int, int> m_jd = mesh.hex_loc_indxs.at(jd);
	//				int j_xi = get<0>(m_jd);
	//				int j_eta = get<1>(m_jd);
	//				int j_zeta = get<2>(m_jd);

	//				d_Nd_xi[jd] = gll::func.dl[order - 1][j_xi][i_xi] * gll::func.l[order - 1][j_eta][i_eta] * gll::func.l[order - 1][j_zeta][i_zeta];
	//				d_Nd_eta[jd] = gll::func.l[order - 1][j_xi][i_xi] * gll::func.dl[order - 1][j_eta][i_eta] * gll::func.l[order - 1][j_zeta][i_zeta];
	//				d_Nd_zeta[jd] = gll::func.l[order - 1][j_xi][i_xi] * gll::func.l[order - 1][j_eta][i_eta] * gll::func.dl[order - 1][j_zeta][i_zeta];
	//			}

	//			// compute J
	//			for (int jd = 0; jd < nodes; jd++)
	//			{

	//				// J[0] dx/d\xi	
	//				// J[1] dy/d\xi
	//				// J[2] dz/d\xi
	//				// J[3] dx/d\eta	
	//				// J[4] dy/d\eta
	//				// J[5] dz/d\eta
	//				// J[6] dx/d\zeta	
	//				// J[7] dy/d\zeta
	//				// J[8] dz/d\zeta

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
	//			assert(abs(detJ) > 1e-14);
	//			Eigen::Matrix3d inv_J = J.inverse();

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
	//			cblas_dgemm(layout, nontrans, nontrans, Ddim, Bcols, Ddim, std::abs(detJ), D.data(), Ddim, B.data(), Ddim, 0, Z.data(), Ddim);
	//			/*	std::ofstream out;
	//				out.open("Z = " + std::to_string(elem_id) + "_" + std::to_string(id) + ".txt");
	//				out << "elem_id = " << elem_id << std::endl << Z << std::endl;
	//				out.close();
	//			*/	cblas_dgemm(layout, trans, nontrans, Bcols, Bcols, Ddim, 1.0, B.data(), Ddim, Z.data(), Ddim, 1, A.data(), Bcols);

	//			/*out.open("A = " + std::to_string(elem_id) + "_" + std::to_string(id) + ".txt");
	//			out << "elem_id = " << elem_id << std::endl << A << std::endl;
	//			out.close();*/
	//		}

	//		std::ofstream out;
	//		out.open("A_" + std::to_string(elem_id) + ".txt");
	//		out << "elem_id = " << elem_id << std::endl << A << std::endl;
	//		out.close();
	//		Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
	//		auto rank = lu_decomp.rank();
	//		std::cout << "elem_id = " << elem_id << " rank = " << rank << std::endl;

	//		for (int id = 0; id < nodes; id++)
	//		{
	//			std::vector<int> glob_index(nodes);
	//			for (int q = rows[idx[id]]; q < rows[idx[id] + 1]; q++)
	//			{
	//				for (int kd = 0; kd < nodes; kd++)
	//				{
	//					if (cols[q] == idx[kd]) glob_index[kd] = q;
	//				}
	//			}
	//			for (int kd = 0; kd < nodes; kd++)
	//			{
	//				//K[blocksize * blocksize * glob_index[kd] + 0] += A(blocksize * id + 0, blocksize * kd + 0);
	//				//K[blocksize * blocksize * glob_index[kd] + 1] += A(blocksize * id + 0, blocksize * kd + 1); //[(1 + 2 * k) * Bcols + 0 + 2 * jd];
	//				//K[blocksize * blocksize * glob_index[kd] + 2] += A(blocksize * id + 0, blocksize * kd + 2);
	//				//															   
	//				//K[blocksize * blocksize * glob_index[kd] + 3] += A(blocksize * id + 1, blocksize * kd + 0);
	//				//K[blocksize * blocksize * glob_index[kd] + 4] += A(blocksize * id + 1, blocksize * kd + 1);
	//				//K[blocksize * blocksize * glob_index[kd] + 5] += A(blocksize * id + 1, blocksize * kd + 2);
	//				//															   
	//				//K[blocksize * blocksize * glob_index[kd] + 6] += A(blocksize * id + 2, blocksize * kd + 0);
	//				//K[blocksize * blocksize * glob_index[kd] + 7] += A(blocksize * id + 2, blocksize * kd + 1);
	//				//K[blocksize * blocksize * glob_index[kd] + 8] += A(blocksize * id + 2, blocksize * kd + 2);

	//				K[blocksize * blocksize * glob_index[kd] + 0] += A(blocksize * kd + 0, blocksize * id + 0);
	//				K[blocksize * blocksize * glob_index[kd] + 1] += A(blocksize * kd + 1, blocksize * id + 0);
	//				K[blocksize * blocksize * glob_index[kd] + 2] += A(blocksize * kd + 2, blocksize * id + 0);

	//				K[blocksize * blocksize * glob_index[kd] + 3] += A(blocksize * kd + 0, blocksize * id + 1);
	//				K[blocksize * blocksize * glob_index[kd] + 4] += A(blocksize * kd + 1, blocksize * id + 1);
	//				K[blocksize * blocksize * glob_index[kd] + 5] += A(blocksize * kd + 2, blocksize * id + 1);

	//				K[blocksize * blocksize * glob_index[kd] + 6] += A(blocksize * kd + 0, blocksize * id + 2);
	//				K[blocksize * blocksize * glob_index[kd] + 7] += A(blocksize * kd + 1, blocksize * id + 2);
	//				K[blocksize * blocksize * glob_index[kd] + 8] += A(blocksize * kd + 2, blocksize * id + 2);

	//			}
	//		}


	//	}

	//	//print_matrix("C:/WD/octave/K.txt", 3, K, rows, cols);
	//}

	void solve(const fc& fcase) 
	{
		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> K;
		buildFullGlobalMatrixStruct(fcase.mesh, rows, cols);
		int breakpoint = 0;
	}
}
