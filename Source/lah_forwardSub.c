#include "lah.h"

#ifdef HAVE_LAPACK /* Use a LAPACK package to do the heavy work */
/*#include <cblas.h>*/
#include <lapacke.h>
lah_Return lah_forwardSub(lah_mat *B, lah_mat const *L)
{
    if ( L == NULL ||  B == NULL || L->nR != B->nR
        || B->nR != L->nC)
    {
        return lahReturnParameterError;
    }
    
    /*
    matrix_layout = LAH_LAPACK_LAYOUT;
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
    
    return (!POTRS(LAH_LAPACK_LAYOUT, 'L', L->nR, B->nC, 
                   L->data, LAH_LEADING_DIM(L), B->data, LAH_LEADING_DIM(B))) ? lahReturnOk : lahReturnExternError;
}
    
#else /* fallback to primitive linear algebra subroutines */
#include <stdlib.h>
#include <math.h>
lah_Return lah_forwardSub(lah_mat *B, const lah_mat  *L)
{
    lah_index i, j, k;
    lah_value *diag, *row, *Aij, *bi, *bj; 

    if ( L == NULL ||  B == NULL || L->nR != B->nR
        || B->nR != L->nC)
    {
        return lahReturnParameterError;
    }

    /* loop over cols of x and b */
    for (k = 0; k < B->nC; k++)
    {
        /* initialize pointers for i-loop*/
        bi = B->data + k * B->incCol;
        diag = L->data;
        row = L->data;
        /* loop over rows of x*/
        for (i = 0; i < B->nR; i++)
        {
            /* initialize pointers for j-loop */
            bj = B->data + k * B->incCol;
            Aij = row;

            /* substitute up to (excluding) current x */
            for (j = 0; j < i; j++)
            {
                *bi -= *Aij * *bj;
                
                /* End of Loop: Increment Pointers */
                bj += B->incRow;
                Aij += L->incCol;
            }
            /* divide xi by diagonal element of L */
            if (fabs(*diag) > 0.00001)
                *bi /= *diag;
            else
                return lahReturnMathError;
            
            /* End of Loop: Increment Pointers */
            diag += L->incRow + L->incCol;
            bi += B->incRow;
            row += L->incRow;
        }
    }
    return lahReturnOk;
}

#endif /* endif primitive implementation */
