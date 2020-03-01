#include "lah.h"

#ifdef HAVE_LAPACK /* Use a LAPACK package to do the heavy work */
#include <lapacke.h>
lah_Return lah_LU(lah_mat *A, lah_index setZero, lapack_int *ipiv)
{
    lah_index i = 0;
    lah_index j = 0;
    lapack_int res = 1;
    if (A == NULL || A->nR != A->nC)
        return lahReturnParameterError;
    
    /*
    matrix_layout = LAH_LAPACK_LAYOUT(A);
    lda = (matrix_layout == LAPACK_COL_MAJOR) ? A->nR : A->nC;
    */
    res = GETRF(LAH_LAPACK_LAYOUT, A->nR, A->nC, A->data, LAH_LEADING_DIM(A), ipiv);
    
    /* Set rest of matrix to zero, if setZero == 1 */
    if (setZero)
    {
        for (i = 0; i < A->nR; i++)
        {
            for(j = i + 1; j < A->nC; ++j)
                LAH_ENTRY(A, i , j) = 0.0;
        }
    }
    return (res == 0) ? lahReturnOk : lahReturnExternError;
}
#endif
 
