#include "lah.h"

/* Adds to matrices sucht that C = alpha * A + beta * B
 * Adding matrices is usually memory band width limited
 * which means there is no performance gain in using BLAS
 * routines. */
lah_Return lah_matAdd(lah_mat *C, 
                      lah_value alpha, lah_mat const *A,
                      lah_value beta, lah_mat const *B)
{
    lah_value *value1 = A->data;
    lah_value *value2 = B->data;
    lah_value *res = C->data;
    lah_index i, j;
    
    if (C == NULL || A == NULL || B == NULL 
        || A->nR != B->nR || A->nC != B->nC  
        || A->nR !=C->nR || A->nC != C->nC)
    {
        return lahReturnParameterError;
    }

    for (i = 0; i < C->nC; i++)        
    {
        value1 = A->data + i * A->incCol;
        value2 = B->data + i * B->incCol; 
        res = C->data + i * C->incCol;
        for (j = 0; j < C->nR; j++)        
        {
            *res = alpha * *value1 + beta * *value2;
            value1 += A->incRow;
            value2 += B->incRow; 
            res += C->incRow;
        }
        
    }
    return lahReturnOk;
}

