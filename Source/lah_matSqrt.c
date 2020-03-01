#include "lah.h"

#ifdef HAVE_LAPACK
#include <lapacke.h>
#include <math.h>

/* Input Matrix P has to be symmetric positive semi-definite, 
 * e.g. real covarianz matrix */
lah_Return lah_matSqrt(lah_mat *P, lah_mat *Psqrt, lah_mat *Work)
{
    lah_Return res = lahReturnOk;
    lah_index i, j;
    lah_value temp;
    lah_value *S;
    lah_mat *Z;
    
    if (P == NULL || P->nR != P->nC
        || Psqrt == NULL || P->nR != Psqrt->nR || P->nC != Psqrt->nC)
    {
        return lahReturnParameterError;
    }
    if (LAH_ISTYPE(P, lahTypeDiagonal))
    {
        for (j = 0; j < P->nC; j++)
            LAH_ENTRY(P, j, j) = sqrt(LAH_ENTRY(P, j, j));
        return lahReturnOk;
    }
    
    /* allocate space for the output parameters and workspace arrays */
    S = malloc(P->nR * sizeof(lah_value));
    Z = lah_matAlloc(P->nR, P->nR, 1);
    res = lah_eigenValue(P, S, Z);
    if(res != lahReturnOk)
        return res;
   
    /* Calculate Workspace = Z*D */
    for (j = 0; j < P->nC; ++j) 
    {
        temp = S[j];
        if (temp < 0)
            return lahReturnMathError; /*Not positive semi definite*/
        else             
            temp = sqrt(temp);
        for (i = 0; i< P->nR; ++i) 
        {
            LAH_ENTRY(Work, i, j) = LAH_ENTRY(Z, i, j) * temp;
        }
    }

    /* Calculate End result P_sqrt = W*Z^T */
    res = lah_matMul(lahNorm, lahTrans, 0, 1, Psqrt, Work, Z);
    
    lah_matFree(Z);
    free(S);
    return res;
}
#endif
