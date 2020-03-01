/*******************************************************************/
/** Test of the LU factorization                                  **/
/*******************************************************************/

#include "lah.h"
#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>

#define PRECISION 4

static lah_mat *example1;
static lah_mat *solution1;
static lah_index *ipiv;

void test_setup(void) 
{
    lah_index nipiv;
    example1 = lah_matLoad("./test_matrices/TestMatrix_LU1.dat"); 
    nipiv = LAH_MIN(example1->nC, example1->nR);
    solution1 = lah_matLoad("./test_matrices/SolMatrix_LU1_U.dat");  
    if (example1 == NULL || solution1 == NULL)
        mu_fail("Could not load matrices\n");
    ipiv = malloc(nipiv * sizeof(lah_index)); 
}

void test_teardown(void) 
{
    lah_matFree(example1);
    lah_matFree(solution1);
    free(ipiv);
}

/*** Test Definitions **********************************************/
/* Example1: Chol
 * ( 0  5  5 )     ( 0 0 1 ) ( 1.0     0.0     0.0 )  ( 6.0 8.0     8.0      )
 * ( 2  9  0 ) --> ( 0 1 0 ) ( 0.33333 1.0     0.0 )  ( 0.0 6.33333 -2.66667 )
 * ( 6  8  8 )     ( 1 0 0 ) ( 0.0     0.78947 1.0 )  ( 0.0 8.0     7.10526  )
 */
MU_TEST(test_LU_example1)
{
    lah_matPrint(example1, 1);
    mu_check(lahReturnOk == lah_LU(example1, 0, ipiv));
    mu_assert_mat_eq(solution1, example1, PRECISION);
    mu_assert_mat_eq(solution1, example1, PRECISION);
}

MU_TEST_SUITE(test_LU_suite)
{   
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	MU_RUN_TEST(test_LU_example1);
}
/*******************************************************************/

int main()
{
    MU_RUN_SUITE(test_LU_suite);
    MU_REPORT();
    
    return 0;
}
