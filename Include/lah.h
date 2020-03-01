#ifndef _LINEAR_ALGEBRA_HELPERS_H_
#define _LINEAR_ALGEBRA_HELPERS_H_
#include <stddef.h>

#ifdef HAVE_LAPACK
#include <lapacke.h>
#include <cblas.h>
#endif

/* Define constants */
#define LAH_PI 3.14159265358979323846

/*--- Data types ------------------------------------------------------------*/

#ifdef HAVE_LAPACK
typedef size_t lah_index;
#else
typedef size_t lah_index;
#endif

/* Choose type of data */
#ifdef USE_SINGLE  /* Single or double precission */
    typedef float lah_value;
#else
    typedef double lah_value;
#endif


/* Automatically choose BLAS/Lapack routines based on data type */
/* If not OpenBLAS (cblas, lapacke) is used, 
 * here is the place to account for another calling convention */
#ifdef HAVE_LAPACK
    #ifdef USE_SINGLE
        #define GEMM cblas_sgemm
        #define POTRF LAPACKE_spotrf
        #define POTRS LAPACKE_spotrs
        #define GETRF LAPACKE_sgetrf
        #define GETRS LAPACKE_sgetrs
        #define GEQRF LAPACKE_sgeqrf
        #define SYEVR LAPACKE_ssyevr
        #define GESDD LAPACKE_sgesdd
        #define LAMCH LAPACKE_slamch
    #else 
        #define GEMM cblas_dgemm
        #define POTRF LAPACKE_dpotrf /* Chol Decomposition */
        #define POTRS LAPACKE_dpotrs /* Solving Triangular System */
        #define GETRF LAPACKE_dgetrf /* LU Decomposition */
        #define GETRS LAPACKE_dgetrs /* Solving Triangular System */
        #define GEQRF LAPACKE_dgeqrf /* QR factorization */
        #define SYEVR LAPACKE_dsyevr /* SVD for matSqrt (symmetric matrix) */
        #define GESDD LAPACKE_dgesdd /* General SVD with divide and conquer */
        #define LAMCH LAPACKE_dlamch /* Helper to obtain machine precision */
    #endif /* end Data_Types*/
#endif /* end HAVE_LAPACK*/

/* Check if called from C++ */
#ifdef __cplusplus
	extern "C" {
#endif

/*--- Definitions for Linear Algebra subroutines*/
/* Specifies whether Row- or ColMajor */
typedef enum
{
    lahRowMajor=101,
    lahColMajor=102
} lah_layout;

#define LAH_LAYOUT lahRowMajor
/*#define LAH_LAYOUT(A) (((A)->matType & lahTypeColMajor) ? lahColMajor : lahRowMajor) */
#ifdef HAVE_LAPACK
    #if LAH_LAYOUT == lahColMajor
        #define LAH_CBLAS_LAYOUT CblasColMajor
        #define LAH_LAPACK_LAYOUT LAPACK_COL_MAJOR
        #define LAH_LEADING_DIM(A) ((A)->nR)
    #else
        #define LAH_CBLAS_LAYOUT CblasRowMajor
        #define LAH_LAPACK_LAYOUT LAPACK_ROW_MAJOR
        #define LAH_LEADING_DIM(A) ((A)->nC)
    #endif
#endif



/* Option for LA operations: Transpose Input before Computation or not? */
typedef enum 
{
    lahTrans = 'T',
    lahNorm = 'N'
} lah_MatOp;

/* Error Types definition*/
typedef enum
{
    lahReturnOk = 0,          /* Success */
    lahReturnParameterError,  /* Wrong Parameters */
    lahReturnMathError,       /* General Math Error, either error while 
                                 * computing or prerequisit not met */
    lahReturnFPError,          /* Error when evaluating function pointer */
    lahReturnExternError          /* Error when evaluating function pointer */
} lah_Return;

typedef enum
{
    /* Properties, which are known in matAlloc */
    lahTypeColMajor = 0x01,     /* 0000 0001 */
    lahTypeNoData = 0x02,       /* 0000 0010 */
    lahTypeVector = 0x04,       /* 0000 0100 */
    lahTypeSquare = 0x08,       /* 0000 1000 */
    
    /* Properties, which are not known in matAlloc */
    lahTypeDiagonal = 0x10,     /* 0001 0000 */
    lahTypeSymmetric = 0x20,    /* 0010 0000 */
    lahTypePositive = 0x40,     /* 0100 0000 */
    lahTypeOrthogonal = 0x80    /* 1000 0000 */
} lah_MatType;

/* lah_mat: Matrix format */
typedef struct 
{
    lah_index nC;          /* Number of Columns */
    lah_index nR;          /* Number of Rows */
    lah_index incRow;      /* increment of rows */
    lah_index incCol;      /* increment of columns */
    lah_value *data;       /* pointer to data */
    lah_MatType matType;   /* flags for several properties of Matrix */
} lah_mat;

/*---------------------------------------------------------------------------*/

/*--- Common Linear Algebra helpers ----------------------------------------*/
/* If provided, these functions are using a CBLAS/LAPACKE library           */
/* else they use a simple handwritten defintion                             */
/*--------------------------------------------------------------------------*/

/* Computes the general Matrix mulitplikation of form:
 * C = alpha*C + beta*op(A)*op(B) 
 * similar to the BLAS routine ?gemm. 
 * op(A) = A for transA=koopNorm or op(A) = At for transA = koopTrans */
lah_Return lah_matMul(lah_MatOp transA, lah_MatOp transB, 
                        lah_value alpha, lah_value beta,
                        lah_mat *C, lah_mat const *A, lah_mat const *B);

/* Adds two matrices C = alpha*A + beta*B */
lah_Return lah_matAdd(lah_mat *C, lah_value alpha, lah_mat const *A, 
                      lah_value beta, lah_mat const *B);

/* This will return a lah_mat struct At with inverted order 
 * (RowMajor <-> ColMajor). If At is not needed anymore, call free(At)      
 * NOT lah_matFree(At), as this will also delete A->data                   */
lah_mat *lah_matTrans(lah_mat const *A);

/* Computes the Cholesky factorization of the form P = L * Lt               
 * where P is a symmetric positive-definite matrix, L is a lower triangular 
 * matrix and Lt denotes its transpose.                                     */
lah_Return lah_chol(lah_mat *P, lah_index setZero);

/* Computes a rank-1 Update to a Cholesky factorization.
 * Let L is the Cholesky factor of a spd matrix P = L * Lt, then
 * the update is the Cholesky factor of the matrix P_new = P + x * x' */
lah_Return lah_cholUpdate(lah_mat *L, lah_mat *x, lah_value alpha);

/* Computes the Forward Substitution of the linear equation system L * X = B
 * where L is a lower triangular matrix, 
 * x and b are matrices with the same number of columns.                    */
lah_Return lah_forwardSub(lah_mat *B, lah_mat const *L);

/* Computes an Update to a matrix C of form D = A*C*A'+B  */
/* workspace needed for special case C = A*C*A'+B         */
lah_Return lah_matUpdate(lah_mat *D, lah_mat const *C, 
                           lah_mat const *A, lah_mat const *B, 
                           lah_mat *W);
/*--------------------------------------------------------------------------*/

/*--- Function only defined for HAVE_LAPACK --------------------------------*/

#ifdef HAVE_LAPACK

/* Computes the LU factorization of the form A = L * U, 
 * as an alternative to Cholesky for general matrices*/
lah_Return lah_LU(lah_mat *A, lah_index setZero, lapack_int *ipiv);

/* Solves the linear equation system L * X = B
 * arising from the LU factorization */
/*lah_Return lah_solveLU(lah_mat *X, const lah_mat *L, 
                       const lah_mat *B, lah_index *ipiv);
*/
lah_Return lah_solveLU(lah_mat *B, const lah_mat *L, lapack_int *ipiv);

/* Computes the QR factorization of the form A = L * U, 
 * as an alternative to Cholesky for general matrices*/
lah_Return lah_QR(lah_mat *A, lah_value *tau);

/* Computes the eigen values and vectors of a symmetric matrix P*/
lah_Return lah_eigenValue(lah_mat *P, lah_value *S, lah_mat *Z);

/* Computes the Sqrt of a positive Matrix by computing the eigen values */
lah_Return lah_matSqrt(lah_mat *P, lah_mat *Psqrt, lah_mat *Work);

/*General Singular Value Decomposition */
lah_Return lah_SVD(lah_mat *A, lah_mat *U, lah_value *S, lah_mat *Vt);
#endif

/*--------------------------------------------------------------------------*/

/*--- Helper functions -----------------------------------------------------*/

/*--- Utility functions ---*/
/* allocation for the Matrix struct lah_mat */
lah_mat* lah_matAlloc(lah_index nC, lah_index nR, lah_index AllocData);

/* Free all data used bei the Matrix struct lah_mat */
lah_mat* lah_matFree(lah_mat *mat);

/* load a matrix from a file */ 
lah_mat* lah_matLoad(const char *fname);

/* Print a lah_mat Matrix struct (used for debugging) */
lah_Return lah_matPrint(lah_mat const *A, lah_index brief);

/*--- Noise Generators ---*/
/* Returns a random number from a Gaussian distribution */
lah_value lah_gaussNoise();

/* Returns a random number from a Lorentz (Cauchy) distribution */
lah_value lah_lorentzNoise();

/*--- Property checks ---*/
/* Checks if matrix is symmetric. 
 * It compares with $(precision) significant digits */
lah_index lah_isSymmetric(lah_mat *A, lah_index precision);

/* Checks if matrix is positive definite. It tries to compute a Cholesky 
 * factorization, which only works for positive definite matrices*/
lah_index lah_isPositive(lah_mat *A);


/*--- Construct different Matrix types ---*/
/* Construct a diagonal matrix from n parameters */
lah_mat *lah_matConstructDiag(lah_index nC, const lah_value *parameter_vec);

/* Construct a symmetric Matrix from n(n+1)/2 parameters */
lah_mat *lah_matConstructSy(lah_index nC, const lah_value *parameter_vec);

/* Construct a orthogonal Matrix from n(n+1)/2 parameters */
lah_mat *lah_matConstructOrtho(lah_index nC, lah_value *parameter_vec);

/* Construct a positive Matrix from n(n+3)/2 parameters */
lah_mat *lah_matConstructPo(lah_index nC, lah_value *parameter_vec);

/* Construct a positive Matrix the simple way from a nxn Maxtrix */
lah_mat *lah_matConstructPo_simple(lah_mat *A);


/*--- Macro definitions ----------------------------------------------------*/

/* Macro to access element of Matrix A[i,j] */
#define LAH_ENTRY(A, i, j) ((A)->data[(i) * (A)->incRow + (j) * (A)->incCol])
#define LAH_MAX(A, B) (((A) > (B)) ? (A) : (B))
#define LAH_MIN(A, B) (((A) < (B)) ? (A) : (B))
#define LAH_SIGN(A) ((A) < 0) ? (-1) : (1))

/*--- Property MACROS---*/
#define LAH_ISTYPE(A, Type) ((A)->matType & (Type))        /* AND */
#define LAH_SETTYPE(A, Type) ((A)->matType |= (Type))      /* OR */
#define LAH_UNSETTYPE(A, Type) ((A)->matType &= (~Type)) /* NAND */
#define LAH_TOGGLETYPE(A, Type) ((A)->matType ^= (Type))   /* XOR */

/*--- LAPACK_ONLY HEADERS ---*/
#ifdef HAVE_LAPACK

#endif 



/*--------------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif

#endif
