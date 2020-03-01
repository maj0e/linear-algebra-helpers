#include "lah.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

lah_index lah_isSymmetric(lah_mat * A, lah_index precision)
{
    lah_index i, j;
    lah_value epsilon = pow(10, (-1) * precision);
    for (j = 0; j < A->nC; j++)
    {    
        for(i = j+1; j < A->nR; j++)
        {
            if (fabs(LAH_ENTRY(A, i, j) - LAH_ENTRY(A, j, i)) > epsilon)
                return 0;
        }
    }
    return 1;
}

lah_index lah_isPositive(lah_mat * A)
{
    lah_mat *Atemp;
    lah_Return result;
    if (!LAH_ISTYPE(A, lahTypeSquare))
        return 0;
    
    Atemp = lah_matAlloc(A->nC, A->nR, 1);
    
    /* Copy data to prevent changes in A */
    memcpy(Atemp->data, A->data, A->nR * A->nC * sizeof(lah_value));
    
    result = lah_chol(Atemp, 0);
    
    lah_matFree(Atemp);
    if (result != lahReturnOk)
        return 0;
    else 
        return 1;
}


