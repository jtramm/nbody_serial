#include "nbody_header.h"

#ifdef MPI
void run_parallel_problem(int nBodies, double dt, int nIters, char * fname)
{
}

void compute_forces_multi_set(Body * local, Body * remote, double dt, int n)
{
}

void parallel_randomizeBodies(Body * bodies, int nBodies_per_rank, int mype, int nprocs)
{
}

// Opens MPI file, and writes header information (nBodies, iterations)
MPI_File initialize_IO(long nBodies, long nIters, char * fname, int mype)
{
	MPI_File fh;
	return fh;
}

// Writes all particle locations for a single timestep
void distributed_write_timestep(Body * local_bodies, long nBodies_per_rank, int timestep, int nIters, int nprocs, int mype, MPI_File * fh)
{
}
#endif
