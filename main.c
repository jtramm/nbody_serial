#include "nbody_header.h"

int main(int argc, char* argv[])
{
	// Input Parameters
	long nBodies = 1000;
	double dt = 0.2; 
	int nIters = 1000;
	int nthreads = 1;
	char * fname = "nbody.dat";

	if( argc != 5 )
	{
		printf("Usage: ./nbody_serial <number of bodies> <number of iterations> <timestep length (dt)> <number of OpenMP threads per rank>\n");
		return 1;
	}

	nBodies = atol(argv[1]);
	nIters = atoi(argv[2]);
	dt = atof(argv[3]);
	nthreads = atoi(argv[4]);

	// Set number of OMP threads if necessary
	#ifdef OPENMP
	omp_set_num_threads(nthreads);
	#endif

	// Initialize MPI
	#ifdef MPI
	MPI_Init(&argc, &argv);
	#endif
	
	// Print Inputs
	print_inputs(nBodies, dt, nIters, nthreads);

	// Run Problem
	#ifdef MPI
	run_parallel_problem(nBodies, dt, nIters, fname);
	#else
	run_serial_problem(nBodies, dt, nIters, fname);
	#endif

	// Finalize MPI
	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
