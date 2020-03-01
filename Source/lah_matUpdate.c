#include "lah.h"

/* Computes an Update to a matrix C of form D = A*C*A'+B  */
lah_Return lah_matUpdate(lah_mat *D, lah_mat const *C, 
                         lah_mat const *A, lah_mat const *B, 
                         lah_mat *W)
{
    lah_Return result;
    result = lah_matMul(lahNorm, lahTrans, 0.0, 1.0, W, C, A);
    result += lah_matMul(lahNorm, lahNorm, 0.0, 1.0, D, A, W);
    result += lah_matAdd(D, 1.0, D, 1.0, B);
    
    return (result == lahReturnOk) ? lahReturnOk : lahReturnMathError;
}

