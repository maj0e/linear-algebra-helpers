/*******************************************************************/
/** Test of the Singular Value decomposition                      **/
/** Only available for HAVE_LAPACK=1                              **/
/*******************************************************************/

#include "lah.h"
#include "minunit.h"
#include <stdio.h>

#define PRECISION -1

lah_mat *example1;
lah_mat *example1_check;
lah_mat *solution1_U;
lah_mat *solution1_Vt;

lah_value *S;
lah_mat *U;
lah_mat *Vt;
lah_mat *work;

void test_setup(void) 
{
    example1 = lah_matLoad("./test_matrices/TestMatrix_SVD1.dat");
    example1_check = lah_matLoad("./test_matrices/TestMatrix_SVD1.dat"); 
    solution1_U = lah_matLoad("./test_matrices/SolMatrix_SVD1_U.dat");
    solution1_Vt = lah_matLoad("./test_matrices/SolMatrix_SVD1_Vt.dat");
    if (example1 == NULL || solution1_U == NULL || solution1_Vt == NULL)
        mu_fail("Could not load matrices\n");
    
    U = lah_matAlloc(solution1_U->nC, solution1_U->nR, 1);
    Vt = lah_matAlloc(solution1_Vt->nC, solution1_Vt->nR, 1);
    S = malloc(solution1_Vt->nR * sizeof(lah_value));
    work = lah_matAlloc(example1->nC, example1->nR, 1);
}

void test_teardown(void) 
{
    lah_matFree(example1);
    lah_matFree(solution1_U);
    lah_matFree(solution1_Vt);
    lah_matFree(example1_check);
}

/*** Test Definitions **********************************************/
/* SVD Example1
 * (  4  12 )  --> ( 0.6   0.8 ) (  20  0 ) (  0.6  0.8 )
 * ( 12  11 )      ( 0.8  -0.6 ) (  0   5 ) ( -0.8  0.6 )
*/
MU_TEST(test_example1)
{
    lah_index i, j;
    mu_check(lahReturnOk == lah_SVD(example1, U, S, Vt));
    mu_assert_mat_eq(solution1_U, U, PRECISION);
    mu_assert_mat_eq(solution1_Vt, Vt, PRECISION);
    
    
    mu_check(fabs(S[0] - 20.0) < 0.001);
    mu_check(fabs(S[1] - 5.0) < 0.001);
    
    /* Sanity check */
    for (i = 0; i < solution1_Vt->nR; i++)
    {  
        for (j = 0; j < solution1_Vt->nC; j++)
        {
            LAH_ENTRY(Vt, i, j) = S[i] * LAH_ENTRY(Vt, i, j);
        }
    }
    
    lah_matMul(lahNorm, lahNorm, 0.0, 1.0, work, U, Vt);
    mu_assert_mat_eq(example1_check, work, PRECISION);
}

MU_TEST_SUITE(test_chol_suite)
{   
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	MU_RUN_TEST(test_example1);
}
/*******************************************************************/

int main()
{
    MU_RUN_SUITE(test_chol_suite);
    MU_REPORT();
    
    return 0;
}
