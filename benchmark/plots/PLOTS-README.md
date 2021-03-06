SLAPS Performance Plots
====

The plots can be generated by simply running the python scripts (with Python 3). Requires matplotlib.

The data was collected on NERSC's Cori Haswell machine, using one node per rank. The CPU speed was limited to 2.3 GHz. The tests were run interactively using the benchmarking scripts in `../slaps` and `../petsc` directories. 

PETSc was built with version 3.9.1, with the options

```--with-scalar-type=complex --with-petsc-arch=haswell-opt CC=cc CXX=CC FC=ftn --with-debugging=0 COPTFLAGS= FOPTFLAGS= CXXOPTFLAGS= --with-shared-libraries=1 LDFLAGS=-dynamic
```

using `PrgEnv-intel/6.0.4` on Cori Haswell.

SLAPS was built as is specified in the benchmark makefile, using `PrgEnv-gnu/6.0.4` and `upcxx/2018.3.0` on Cori Haswell.