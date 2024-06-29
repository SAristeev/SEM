#include "algebra.h"

#include <iostream>

#include <vector>

#define AMGCL_NO_BOOST

#include <amgcl/preconditioner/schur_pressure_correction.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/make_block_solver.hpp>
#include <amgcl/preconditioner/schur_pressure_correction.hpp>
#include <amgcl/value_type/eigen.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/spai1.hpp>
#include <amgcl/relaxation/ilu0.hpp>

#include <amgcl/relaxation/as_preconditioner.hpp>


#include <amgcl/solver/preonly.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>


#include <amgcl/profiler.hpp>
#include <amgcl/io/mm.hpp>

#include <mkl.h>

namespace algebra 
{
	void bsr2csr(const int& blocksize, const std::vector<double>& bsr_A, const std::vector<int>& bsr_rows, const std::vector<int>& bsr_cols,
		std::vector<double>& csr_A, std::vector<int>& csr_rows, std::vector<int>& csr_cols) 
	{
		csr_A.resize(bsr_A.size());
		csr_rows.resize(blocksize * (bsr_rows.size() - 1) + 1);
		csr_cols.resize(blocksize * blocksize * bsr_cols.size());

		for (int row = 0; row < bsr_rows.size() - 1; row++)
		{
			int csr_row = blocksize * (bsr_rows[row + 1] - bsr_rows[row]);
			for (int bl = 0; bl < blocksize; bl++)
			{
				csr_rows[3 * row + bl + 1] = csr_rows[3 * row + bl] + csr_row;
			}
		}
		
		int block_row = 0;
		for (int block = 0; block < bsr_cols.size(); block++)
		{
			int block_col = blocksize * bsr_cols[block];
			
			if (block >= bsr_rows[block_row + 1])
			{
				block_row++;
			}
			int bsr_shift = block - bsr_rows[block_row];
			for (int loc_row = 0; loc_row < blocksize; loc_row++)
			{
				int csr_row_shift = csr_rows[blocksize * block_row + loc_row];
				for (int loc_col = 0; loc_col < blocksize; loc_col++)
				{
					int tmp = 0;
					csr_cols[csr_row_shift + blocksize * bsr_shift + loc_col] = block_col + loc_col;
					csr_A[csr_row_shift + blocksize * bsr_shift + loc_col] = bsr_A[blocksize * blocksize * block + blocksize * loc_row + loc_col];
				}
			}
			
		}
	}

	void solve_pardiso(const int& blocksize, const std::vector<double>& A, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x) 
	{
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

	void solve_amgcl(const int& blocksize, const std::vector<double>& A, const std::vector<int>& rows, const std::vector<int>& cols, const int& nrhs, const std::vector<double>& b, std::vector<double>& x) 
	{
		assert(blocksize == 1);
		int n = rows.size() - 1;
		auto Matrix = std::tie(n, rows, cols, A);

		typedef amgcl::backend::builtin<double> Backend;
		typedef amgcl::make_solver<
			amgcl::amg<
			Backend,
			amgcl::coarsening::ruge_stuben,
			amgcl::relaxation::gauss_seidel
			>,
			amgcl::solver::cg<Backend>
		> Solver;
		
		Solver::params prm;
		//prm.solver.M = 100;
		prm.precond.npre = 4;
		prm.precond.npost = 4;
		//prm.precond.
		//prm.precond.ncycle = 4;

		prm.solver.tol = 1e-9;
		prm.solver.maxiter = 10000;
		prm.solver.verbose = true;

		// The profiler:
		amgcl::profiler<> prof("SEM");
		// Initialize the solver with the system matrix:
		prof.tic("setup");
		Solver solve(Matrix, prm);
		prof.toc("setup");

		// Show the mini-report on the constructed solver:
		std::cout << solve << std::endl;


		// Solve the system with the zero initial approximation:
		int iters;
		double error;

		prof.tic("solve");
		std::tie(iters, error) = solve(Matrix, b, x);
		prof.toc("solve");

		// Output the number of iterations, the relative error,
		// and the profiling data:
		std::cout << "Iters: " << iters << std::endl
			<< "Error: " << error << std::endl
			<< prof << std::endl;

	}
}