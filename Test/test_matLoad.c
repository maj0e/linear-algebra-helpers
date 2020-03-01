#include "lah.h"
#include "minunit.h"
#include <stdlib.h>
#include <stdio.h>

/*--- Define global data for tests ---*/
lah_mat *A_correct;
lah_mat *A = NULL;

lah_value Adata_rmaj[] = { 1,  2,  3,  4,
                           5,  6,  7,  8,
                           9, 10, 11, 12,
                          13, 14, 15, 16,};

lah_value Adata_cmaj[] = { 1,  5,  9, 13,
                         2,  6, 10, 14,
                         3,  7, 11, 15,
                         4,  8, 12, 16,};

/*--- Define setup and finish of test suite ---*/

void test_setup(void) 
{
    if (LAH_LAYOUT == lahColMajor)
    {    
        A_correct= lah_matAlloc(4 ,4, 0);
        A_correct->data = Adata_cmaj;
    }
    else
    {
        A_correct = lah_matAlloc(4 ,4, 0);
        A_correct->data = Adata_rmaj; 
    }
}

void test_teardown(void) 
{
    free(A_correct);
    lah_matFree(A);
}

/*--- Define Tests ---*/
MU_TEST(load_square_matrix) {
    A = lah_matLoad("./test_matrices/TestMatrix_matLoad.dat");
	mu_check(A != NULL);
    mu_assert_mat_eq( A_correct, A, 4) ;
}

/*--- Define Test suite ---*/
MU_TEST_SUITE(test_suite) {
	MU_SUITE_CONFIGURE(&test_setup, &test_teardown);
	MU_RUN_TEST(load_square_matrix);
}

int main()
{
    MU_RUN_SUITE(test_suite);
	
    /* MU_REPORT(); */
	
    return 0;
}
