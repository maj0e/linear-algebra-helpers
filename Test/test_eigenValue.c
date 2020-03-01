/*******************************************************************/
/** Test of the Eigen Decomposition                               **/
/** Only available for HAVE_LAPACK=1                              **/
/*******************************************************************/

#include "lah.h"
#include "minunit.h"
#include <stdio.h>

#define PRECISION -1

lah_mat *example1;
lah_mat *solution1_z;

lah_value *S;
lah_mat *Z;
lah_mat *work;

void test_setup(void) 
{
    example1 = lah_matLoad("./test_matrices/TestMatrix_SVD1.dat");    
    solution1_z = lah_matLoad("./test_matrices/SolMatrix_SVD1_U.dat");
    if (example1 == NULL || solution1_z == NULL)
        mu_fail("Could not load matrices\n");
    
    Z = lah_matAlloc(solution1_z->nC, solution1_z->nR, 1);
    S = malloc(example1->nR * sizeof(lah_value));
}

void test_teardown(void) 
{
    lah_matFree(example1);
    lah_matFree(solution1_z);
    lah_matFree(Z);
    free(S);
}

/*** Test Definitions **********************************************/
/* SVD Example1
 * (  4  12 )  --> ( 0.6   0.8 ) (  20  0 ) (  0.6  0.8 )
 * ( 12  11 )      ( 0.8  -0.6 ) (  0   5 ) ( -0.8  0.6 )
*/
MU_TEST(test_example1)
{
    /* lah_index i, j; */
    mu_check(!lah_eigenValue(example1, S, Z));
    mu_assert_mat_eq(solution1_z, Z, PRECISION);
    
    mu_check(S[0] == 20.0);
    mu_check(S[1] == 5.0);
    
    /* Sanity check */
    /*
    for (i = 0; i < solution1_Vt->nR; i++)
    {  
        for (j = 0; i < solution1_Vt->nC; i++)
        {
            LAH_ENTRY(Vt, i, j) = S[j] * LAH_ENTRY(Vt, i, j);
        }
    }
    
    lah_matMul(lahNorm, lahNorm, 0.0, 1.0, work, U, Vt);
    mu_assert_mat_eq(example1, work, PRECISION);
    */
}

MU_TEST_SUITE(test_EigenValue_suite)
{   
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	MU_RUN_TEST(test_example1);
}
/*******************************************************************/

int main()
{
    MU_RUN_SUITE(test_EigenValue_suite);
    MU_REPORT();
    
    return 0;
}
