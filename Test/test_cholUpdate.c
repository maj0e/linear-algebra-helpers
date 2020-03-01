/*******************************************************************/
/** Test of the Cholesky factorization                            **/
/** Should be tested with all combinations of HAVE_LAPACK,        **/
/** Row-/ColMajor and USE_SINGLE                                  **/
/*******************************************************************/

#include "lah.h"
#include "minunit.h"
#include <stdio.h>

#define PRECISION -1

lah_mat *example1;
lah_mat *example1_comp;
lah_mat * update_vec;

#define NAlPHA 4
lah_value alpha_values[NAlPHA] = {1.0, -0.3, 2.0, 0.0};
lah_value alpha;
char alpha_message[1024];

void test_setup(void) 
{
    example1 = lah_matLoad("./test_matrices/TestMatrix_Chol1.dat");
    example1_comp = lah_matLoad("./test_matrices/TestMatrix_Chol1.dat");
    update_vec = lah_matLoad("./test_matrices/update_vec_cholupdate1.dat"); 
    if (!example1 || !example1_comp || !update_vec)
        mu_fail("Could not load matrices\n");
}

void test_teardown(void) 
{
    lah_matFree(example1);
    lah_matFree(example1_comp);
    lah_matFree(update_vec);
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
    mu_check(lahReturnOk == lah_matMul(lahNorm, lahTrans, 1.0, alpha, 
                                       example1_comp, update_vec, update_vec));
    mu_check(lahReturnOk == lah_chol(example1_comp, 1));
    mu_check(lahReturnOk == lah_chol(example1, 1));
    mu_assert(lahReturnOk == lah_cholUpdate(example1, update_vec, alpha), alpha_message);
    mu_assert_mat_eq(example1_comp, example1, PRECISION);
}


MU_TEST_SUITE(test_chol_suite)
{   
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown);
    for (lah_index i = 0; i < NAlPHA; i++ )
    {	
        alpha = alpha_values[i];
        sprintf(alpha_message, "Error occured for alpha = %g", alpha);
        MU_RUN_TEST(test_example1);
    }
}
/*******************************************************************/

int main()
{
    MU_RUN_SUITE(test_chol_suite);
    MU_REPORT();
    
    return 0;
}
