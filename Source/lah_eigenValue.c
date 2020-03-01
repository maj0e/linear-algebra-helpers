#include "lah.h"

#ifdef HAVE_LAPACK
#include <lapacke.h>

lah_Return lah_eigenValue(lah_mat *P, lah_value *S, lah_mat *Z)
{
    lapack_int res, nEig;
    lapack_int *isuppz;
    if (P == NULL || S == NULL || Z == NULL 
        || P->nR != P->nC)
    {
        return lahReturnParameterError;
    }
    /* allocate space for the output parameters and workspace arrays */
    isuppz = malloc(2 * P->nR * sizeof(lapack_int));

    /* Data allocation */
    /* Sollte A-> data nur dreiecksmatrix enthalten ?*/
    res = SYEVR(LAH_LAPACK_LAYOUT, 'V', 'A', 'L', P->nR, P->data, P->nR, 
                0, 0, 0, 0, LAMCH('S'), &nEig, S, Z->data, P->nR, isuppz);
    
    /* free support array */
    free(isuppz);
    if (res > 0)
        return lahReturnExternError;
    else if (res < 0)
        return lahReturnParameterError;
    else 
        return lahReturnOk;
}    
#endif
