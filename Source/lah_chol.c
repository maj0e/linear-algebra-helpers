#include "lah.h"

#ifdef HAVE_LAPACK /* Use a LAPACK package to do the heavy work */
#include <lapacke.h>
lah_Return lah_chol(lah_mat *P, lah_index setZero)
{
    lah_index i, j;
    lapack_int res;
    if (P == NULL || !LAH_ISTYPE(P, lahTypeSquare))
        return lahReturnParameterError;
    
    res = POTRF(LAH_LAPACK_LAYOUT, 'L' , P->nR, P->data, P->nR);
    
    /* Set rest of matrix to zero, if setZero == 1 */
    if (setZero)
    {
        for (i = 0; i < P->nR; i++)
        {
            for(j = i + 1; j < P->nC; ++j)
            {
                LAH_ENTRY(P, i, j) = 0.0;
            }
        }
    }
    return (res == 0) ? lahReturnOk : lahReturnExternError;
}

#else /* fallback to primitive linear algebra subroutines */
#include <math.h>

lah_Return lah_chol(lah_mat *P, lah_index setZero)
{
    lah_index i, j, k, N;
    lah_value *diag, *value;

    if (P == NULL)
    {
        return lahReturnParameterError;
    }

    N = P->nR;

    /* see https://de.wikipedia.org/wiki/Cholesky-Zerlegung#Pseudocode */
    for (i = 0; i < N; i++)
    {
        /* Pointer to diagonal element */
        diag = P->data + i * (P->incRow + P->incCol);
        /* check if positive */
        if (*diag < 0.00001 )
        {
            return lahReturnMathError;
        }
        /* square root of diagonal element in place */
        *diag = sqrt(*diag);
        
        /* divide elements below diagonal by diagonal element */
        for (k = 1; k < N - i; k++)
        {
            *(diag + k * P->incRow) /= *diag;
        }
        /* compute upper right submatrix */
        for (j = i + 1; j < N; j++)
        {
            value = P->data + j * P->incRow + j * P->incCol;
            for (k = j; k < N; k++)
            {
                *value -= LAH_ENTRY(P, k, i) * LAH_ENTRY(P, j, i);
                value += P->incRow;
            }
        }
    }
    /* Set rest of matrix to zero, if setZero == 1 */
    if (setZero)
    {
        for (i = 0; i < N; i++)
        {
            for(j = i + 1; j < N; ++j)
            {
                LAH_ENTRY(P, i, j) = 0.0;
            }
        }
    }

    return lahReturnOk;
}
#endif /* endif primitive implementation */
