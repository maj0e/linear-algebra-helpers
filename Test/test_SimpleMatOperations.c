#include "lah.h"
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************/
/** Test of the simple Math operations                            **/
/** The LAPACK versions schould yield the same results            **/
/*******************************************************************/

int main()
{
    /* Test should be performed for RowMajor and ColMajor order */
    lah_mat *Bt;
    lah_mat *BtView;
    lah_value Adata_row[6] =  { 1.0,  2.0,  3.0, 
                                 4.0,  5.0,  6.0 };
    lah_value Bdata_row[6] =  { 6.0,  5.0,  4.0,
                                 3.0,  2.0,  1.0 };
    lah_value Adata_col[6] =  { 1.0,  4.0,  
                                 2.0,  5.0,
                                 3.0,  6.0 };
    lah_value Bdata_col[6] =  { 6.0,  3.0,  
                                 5.0,  2.0,
                                 4.0,  1.0 };

               
    /* Create koop_mat structs */
    lah_mat *E = lah_matAlloc(2, 2, 1);
    lah_mat *D = lah_matAlloc(3, 3, 1);
    lah_mat *C = lah_matAlloc(3, 2, 1);
    lah_mat *A = lah_matAlloc(3, 2, 0);
    lah_mat *B = lah_matAlloc(3, 2, 0);
    
    if (LAH_LAYOUT == lahRowMajor)
    {
        A->data = Adata_row;
        B->data = Bdata_row;
    }
    else
    {
        A->data = Adata_col;
        B->data = Bdata_col;
    }

    
    printf("Simple Mat operations begin...\n");
    /* Show Input Matrices */
    printf("Matrix  A: \n");
    lah_matPrint(A,1);

    printf("Matrix  B: \n");
    lah_matPrint(B,1);
    
    /* Test koop_matAdd */
    lah_matAdd(C, 1.0, A, 1.0, B);
    printf("Result of A + B:\n");
    lah_matPrint(C,1); 
    
    lah_matAdd(C, 0.5, A, -2.0, B);
    printf("Result of 0.5 * A - 2 * B:\n");
    lah_matPrint(C,1);
    
    /* Test transpose and transpose view */
    Bt = lah_matAlloc(2, 3, 1);
    /*lah_matTrans(Bt, B); */
    BtView = lah_matTrans(B);
    /*
    printf("Transposed of B:\n");
    lah_matPrint(Bt, 1);
    */
    printf("Transposed View of B (should be the same):\n");
    lah_matPrint(BtView, 1);
    
    
    /* Test Matrix Multiplikation koop_matMul */
    lah_matMul(lahNorm, lahNorm, 0.0, 0.5, E, A, BtView);
    printf("Result of 0.5 * (A * B')\n");
    lah_matPrint(E, 1);

    lah_matMul(lahNorm, lahNorm, 0.0, 0.5, D, BtView, A);
    printf("Result of 0.5 * (B' * A)\n");
    lah_matPrint(D, 1);

    lah_matMul(lahNorm, lahTrans, 0.0, 0.5, E, A, B);
    printf("Result of 0.5 * (A * B') (with Transpose from matMul)\n");
    lah_matPrint(E, 1);

    free(A); /* Don't matFree. Data resides in array */
    free(B); /* Don't matFree. Data resides in array */
    
    free(BtView);
    
    lah_matFree(Bt);
    lah_matFree(C);
    lah_matFree(D);
    lah_matFree(E);
    return 0;


}
