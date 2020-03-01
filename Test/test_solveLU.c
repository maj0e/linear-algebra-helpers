/*******************************************************************/
/** Test of the Forward Substitution                              **/
/** Should be tested with all combinations of HAVE_LAPACK,        **/
/** Row-/ColMajor and USE_SINGLE                                  **/
/** After Testing: ForwardSub(X, A, B ) also works with X==B      **/
/*******************************************************************/
/* Example1: Solve the system of linear equations:
 * (  3   0   0  0 ) (x_1)   (5)     (x_1)   (  1.66667 )
 * ( -1   1   0  0 ) (x_2) = (6) --> (x_2) = (  7.66667 )
 * (  3  -2  -1  0 ) (x_3)   (4)     (x_3)   (-14.3333  )
 * (  1  -2   6  2 ) (x_4)   (2)     (x_4)   ( 50.8333  )
*/
#include "lah.h"
#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>

/*Significant digits to compare floats*/
#define PRECISION 4

lah_mat *A_ex1;
lah_mat *B_ex1;
lah_mat *Sol_ex1;
lah_mat *Res_ex1;
lah_mat *temp;
lah_index *ipiv;
void test_setup(void) 
{
    lah_index nipiv;
    A_ex1 = lah_matLoad("./test_matrices/TestMatrix_ForwardSub1.dat");    
    B_ex1 = lah_matLoad("./test_matrices/TestMatrix_ForwardSub1_rhs.dat");
    Res_ex1 = lah_matLoad("./test_matrices/TestMatrix_ForwardSub1_rhs.dat");
    Sol_ex1 = lah_matLoad("./test_matrices/SolMatrix_ForwardSub1.dat");
    if (A_ex1 == NULL || B_ex1 == NULL || Sol_ex1 == NULL)
        mu_fail("Could not load matrices\n");
    
    temp = lah_matAlloc(Sol_ex1->nC, Sol_ex1->nR, 1);
    nipiv = LAH_MAX(A_ex1->nC, A_ex1->nR);
    ipiv = malloc(nipiv * sizeof(index));
}

void test_teardown(void) 
{
    lah_matFree(A_ex1);
    lah_matFree(B_ex1);
    lah_matFree(Sol_ex1);
    lah_matFree(Res_ex1);
    lah_matFree(temp);
    free(ipiv);
}

MU_TEST(test_example1)
{
    mu_check(!lah_solveLU(Res_ex1, A_ex1, ipiv));
    mu_assert_mat_eq(Sol_ex1, Res_ex1, PRECISION);
    
    /* Double check for consistency*/
    lah_matMul(lahNorm, lahNorm, 0.0, 1.0, temp, A_ex1, Res_ex1);
    mu_assert_mat_eq(B_ex1, temp, PRECISION); 
}

MU_TEST_SUITE(test_fs_suite)
{   
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	MU_RUN_TEST(test_example1);
}

int main()
{
    MU_RUN_SUITE(test_fs_suite);
    MU_REPORT();
    
    return 0;
}
