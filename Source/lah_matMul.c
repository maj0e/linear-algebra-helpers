#include "lah.h"

#ifdef HAVE_LAPACK /* Use a LAPACK package to do the heavy work */
#include <cblas.h>

lah_Return lah_matMul(lah_MatOp transA, lah_MatOp transB,
                      lah_value alpha, lah_value beta, 
                      lah_mat *C, lah_mat const *A, lah_mat const *B)
{
    
    CBLAS_TRANSPOSE transA_cblas = CblasNoTrans;
    CBLAS_TRANSPOSE transB_cblas = CblasNoTrans;
        
    if (C == NULL || A == NULL || B == NULL )
        return lahReturnParameterError;
    
    /*
    order = LAH_CBLAS_LAYOUT(C);
        
    if (order == CblasColMajor)
    {
        lda = A->nR;
        ldb = B->nR;
        ldc = C->nR;
    }
    else 
    { 
        lda = A->nC;
        ldb = B->nC;
        ldc = C->nC;
    }
    */
        
    if (transA != lahNorm)
        transA_cblas = CblasTrans;
    if (transB != lahNorm)
        transB_cblas = CblasTrans;
    
    GEMM(LAH_CBLAS_LAYOUT, transA_cblas, transB_cblas, C->nR, C->nC, A->nC, beta,
         A->data, LAH_LEADING_DIM(A), B->data, LAH_LEADING_DIM(B), alpha, C->data, LAH_LEADING_DIM(C));
        
    return lahReturnOk;
}
#else /* fallback to primitive linear algebra subroutines */
#include <stdlib.h>
/* Computes the general Matrix mulitplikation of form:
 * C = alpha*C + beta*op(A)*op(B) 
 * similar to the BLAS routine ?gemm. 
 * op(A) = A for transA=lahNorm or op(A) = At for transA = lahTrans */
lah_Return lah_matMul(lah_MatOp transA, lah_MatOp transB, 
                      lah_value alpha, lah_value beta,
                      lah_mat *C, const lah_mat *A, const lah_mat *B)
{
    lah_index c, r, i;
    lah_value *res, *value1, *value2;
    lah_value temp;

    if (C == NULL || A == NULL || B == NULL)    
    {
        return lahReturnParameterError;
    }

    /* NOTE: lah_matMul calls itself with an transposed view of 
     * the corresponding kop_mat and sets transA resp. transB 
     * to lahNorm (else we would have an endless recursion) 
    */
    if (transA == lahTrans)
    {
        /* Call lah_matMul again with transposed */
        lah_mat *At = lah_matTrans(A);
        lah_Return result = lah_matMul(lahNorm, transB, alpha,
                                         beta, C, At, B);
        free(At);
        return result;
    }
    if (transB == lahTrans)
    {
        /* Call lah_matMul again with transposed of B */
        lah_mat *Bt = lah_matTrans(B);
        lah_Return result = lah_matMul(transA, lahNorm, alpha,
                                       beta, C, A, Bt);
        free(Bt);
        return result;
    }
    
    if (C->nR != A->nR || C->nC != B->nC || A->nC != B->nR)
    {
        return lahReturnParameterError;
    }

    for (c = 0; c < C->nC; c++)
    {
        for (r = 0; r < C->nR; r++)
        {
            res = C->data + C->incRow * r + c * C->incCol;
            value1 = A->data + r * A->incRow;
            value2 = B->data + c * B->incCol;
            temp = 0;
            for (i = 0; i < A->nC; i++, value1 += A->incCol, value2 += B->incRow)
            {
                temp += *value1 * *value2;
            }
            *res = alpha * *res + beta * temp;
        }
    }

    return lahReturnOk;
}
#endif /* endif primitive implementation */
