/*******************************************************************/
/** Test of the Cholesky factorization                            **/
/** Should be tested with all combinations of HAVE_LAPACK,        **/
/** Row-/ColMajor and USE_SINGLE                                  **/
/*******************************************************************/
/* Example to compute
 * 18  22   54   42           4.24264    0.00000    0.00000    0.00000
 * 22  70   86   62   -->     5.18545    6.56591    0.00000    0.00000
 * 54  86  174  134          12.72792    3.04604    1.64974    0.00000
 * 42  62  134  106           9.89949    1.62455    1.84971    1.39262
*/
#include "lah.h"
#include <stdio.h>
#include <stdlib.h>

int main()
{
    /* A is symmetric: No extra array for ColMajor*/
    lah_value Adata[4*4] =  { 18,  22,  54,  42,
                               22,  70,  86,  62,
                               54,  86, 174, 134,
                               42,  62, 134, 106 };

    lah_mat *A = lah_matAlloc(4, 4, 0);
    A->data = Adata;
        
    printf("Example Matrix A: \n");
    lah_matPrint(A,1);

    /* Do factorization */
    lah_chol(A, 1);

    printf("Factorization of A: \n");
    printf("lower triangle, upper triangle left unchanged) \n");
    lah_matPrint(A,1);

    free(A);
    return 0;
}
