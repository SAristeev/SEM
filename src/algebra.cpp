#include "algebra.h"

#include <vector>

#include <mkl.h>

namespace algebra 
{
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
}