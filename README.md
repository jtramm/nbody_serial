# nbody_serial
A serial solver for the N-Body gravitational force problem in 3D.

This code was written as a starting point for a problem set in
MPCS 51087 at the University of Chicago.

## Compilation
Compilation settings can be configured at the top of the included makefile.
Settings include turning on/off MPI and OpenMP.

$> make

## Runtime Settings

- The total number of bodies to simulate can be set at runtime using the first command line argument of the program.

- The total number of timesteps to run can be set at runtime using the second command line argument.

- The timestep length (delta t) can be set at runtime using the third command line argument.

- The number of threads to run per MPI rank can be set using the fourth command line argument

## Running
example:

```bash
$> make
$> ./nbody 1000 100 0.2 4
```

Will simulate 1000 bodies for 100 timesteps each of length 0.2, using 4 OpenMP threads per rank.

**NOTE:** This code has not yet been parallelized, so using more than one OpenMP thread will not have any effect.

## Parallel Framework

A parallel framework has been included for you in the "parallel.c" file that includes a few function prototypes and a driver function. Feel free to use this to begin your implementation --or-- toss it out an go your own way! You are not required to meet these interfaces, you are free to hack it up however you feel is best.

## Plotting
The included plotter.py python script can be used to generate an animation using the results of the program stored in the ```nbody.dat``` file. The script does not take any arguments. However, when plotting big endian files (like those generated by BG/Q Vesta), you will want to set the "BigEndian" variable to "True" at the top of the script.

You may also want to adjust the view bounds found later in the script depending on how you set up your nbody initial conditions.

On some systems, you may also need to adjust the bitrate, resolution, and DPI settings at the bottom of the script.
