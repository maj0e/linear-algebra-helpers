#include "lah.h"

#ifdef HAVE_LAPACK
#include <lapacke.h>

/* General Singular Value Decomposition*/
/* For now get all eigen values jobz ='A', 
 * Can use this parameter later for truncated version */
lah_Return lah_SVD(lah_mat *A, lah_mat *U, lah_value *S, lah_mat *Vt)
{
    lapack_int res;
    
    if (A == NULL || U == NULL || S == NULL || Vt == NULL 
        || A->nR != A->nC || A->nR != U->nR || A->nC != Vt->nC )
    {
        return lahReturnParameterError;
    }
    /*
    matrix_layout = LAH_LAPACK_LAYOUT(A);
    if (matrix_layout == LAPACK_COL_MAJOR)
    {
        lda = A->nR;
        ldu = U->nR;
        ldv = Vt->nR;
    }
    else 
    { 
        lda = A->nC;
        ldu = U->nC;
        ldv = Vt->nC;
    }
    */

    /* Do the actual SVD: */ 
    res = GESDD(LAH_LAPACK_LAYOUT, 'A', A->nC, A->nR, A->data, LAH_LEADING_DIM(A), S, U->data, LAH_LEADING_DIM(U), Vt->data, LAH_LEADING_DIM(Vt));
    return (res == 0) ? lahReturnOk : lahReturnExternError;
}
#endif
