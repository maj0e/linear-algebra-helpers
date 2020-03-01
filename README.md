# LAH - Linear Algebra Helpers

A small linear algebra library written in C. 

It is closely related to the [kalman filter project](https://gitlab.com/optimal-control/kalman-filter)
and serves two main purposes:

	- It is a wrapper for LAPACK for all functions related to the unscented Kalman filter
	- It implements all functions needed for the (extended) Kalman filter, without LAPACK 

Beside that, it is hopefully useful as a showcase for the usage of BLAS/LAPACK routines, since
information for newcomers can be hard to find.

## Functionality
### General
LAH can be used with column or row major storage format, for compatibility with other libraries (e.g. armadillo).
Therefore the general Matrix looks like this:

```c
/* lah_mat: Matrix format */
typedef struct
{
        lah_index nC;          /* Number of Columns */
        lah_index nR;          /* Number of Rows */
        lah_index incRow;      /* increment of rows */
        lah_index incCol;      /* increment of columns */
        lah_MatType matType;   /* flags for several properties of Matrix */
        lah_value *data;       /* pointer to data */
} lah_mat;                                                       
```
The additional flag `matType` can be used to specify certain properties 
of the matrix (eg. square, symmetric, positive definite etc.),
which can then be used to choose faster LAPACK subroutines or other shortcuts.
Some basic properties (layout, square etc) are set, when the matrix is allocated via `lah_matAlloc`.
Others are set when explicitely constructing such kind of matrices (e.g. `lah_constructPositive`) 
or have to be set manually with ` LAH_SETTYPE(A, Type)` by the user. 


## Getting started
Prerequisites:

	- C compiler (tested with gcc and clang)
	- make
	
Compile static release version with

```
	make OPT=1
```

For the whole functionality openBLAS is needed:

```
	make OPT=1 HAVE_LAPACK=1
```

It's probably possible to use another LAPACK package, 
by setting the appropriate flags in the Makefile and 
change the LAPACK macros in lah.h to account for another calling convention.
