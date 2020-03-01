#include "lah.h"

#ifdef HAVE_LAPACK
#include <lapacke.h>
lah_Return lah_solveLU(lah_mat *B, lah_mat const *L, lapack_int *ipiv)
{
    if ( L == NULL ||  B == NULL || L->nR != B->nR
        || B->nR != L->nC)
    {
        return lahReturnParameterError;
    }
    /*
    if (matrix_layout == LAPACK_COL_MAJOR) 
    {
        lda = L->nR;
        ldb = B->nR;
    }
    else 
    { 
        lda = L->nC;
        ldb = B->nC;
    }
    */
    
    return (!GETRS(LAH_LAPACK_LAYOUT, 'N', B->nR, B->nC, L->data , LAH_LEADING_DIM(L), ipiv, B->data, LAH_LEADING_DIM(B))) ? lahReturnOk : lahReturnExternError;
}
#endif

