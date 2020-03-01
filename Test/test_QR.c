 
/*******************************************************************/
/** Test of the Cholesky factorization                            **/
/** Should be tested with all combinations of HAVE_LAPACK,        **/
/** Row-/ColMajor and USE_SINGLE                                  **/
/*******************************************************************/

#include "lah.h"
#include "minunit.h"
#include <stdio.h>

#define PRECISION 5

lah_mat *example1;
lah_mat *solution1;
lah_value *tau;

lah_order order;

void test_setup(void) 
{
    example1 = lah_matLoad("./test_matrices/TestMatrix_QR1.dat", order);    
    solution1 = lah_matLoad("./test_matrices/SolMatrix_Chol1.dat", order);  
    if (!example1 || !solution1)
        mu_fail("Could not load matrices\n");
    tau = malloc(LAH_MIN(example1->nC, example1-nR) * sizeof(lah_value))
}

void test_teardown(void) 
{
    lah_matFree(example1);
    lah_matFree(solution1);
    free(tau);
}

/*** Test Definitions **********************************************/
/* Example1: Chol
 * 18  22   54   42           4.24264    0.00000    0.00000    0.00000
 * 22  70   86   62   -->     5.18545    6.56591    0.00000    0.00000
 * 54  86  174  134          12.72792    3.04604    1.64974    0.00000
 * 42  62  134  106           9.89949    1.62455    1.84971    1.39262
*/
MU_TEST(test_example1)
{
    mu_check(lahReturnOk == lah_QR(example1, tau));
    lah_
    matPrint(example1, 1);
    mu_assert_mat_eq(solution1, example1, PRECISION);
}

MU_TEST_SUITE(test_chol_suite)
{   
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	MU_RUN_TEST(test_example1);
}
/*******************************************************************/

int main()
{
    order = lahRowMajor;
    printf("Run Cholesky tests for Row Major\n");
    MU_RUN_SUITE(test_chol_suite);
    order = lahColMajor;
    printf("Run Cholesky tests for Col Major\n");
    MU_RUN_SUITE(test_chol_suite);
    MU_REPORT();
    
    return 0;
