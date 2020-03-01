#include "lah.h"
#include <math.h>
lah_Return lah_cholUpdate(lah_mat *L, lah_mat *x, lah_value alpha)
{
    if (L == NULL || x == NULL 
        || L->nR != L->nC || L->nR != x->nR) 
    {
        return lahReturnParameterError;
    }
    
    lah_index i, j;
    lah_value r, c, s, temp;
    for (i = 0; i < x->nR; i++)
    {
        temp = LAH_ENTRY(L, i, i);
        r = sqrt(temp * temp + alpha * x->data[i] * x->data[i]);
        c = r / temp;
        s = x->data[i] / temp;
        LAH_ENTRY(L, i, i) = r;
        
        if (i < x->nR - 1)
        {
            for (j = i + 1; j < x->nR; j ++)
            {
                LAH_ENTRY(L, j, i) = ( LAH_ENTRY(L, j, i) + alpha *  s * x->data[j]) / c;
                x->data[j] = c * x->data[j] - s * LAH_ENTRY(L, j, i);
            }
        }
    }
    return lahReturnOk;
}
