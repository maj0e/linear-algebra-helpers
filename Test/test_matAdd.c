#include "lah.h"
#include "minunit.h"

lah_mat *C_exp;
lah_mat *C_res;
lah_mat *A;
lah_mat *B;
lah_mat *D;

void test_setup(void) 
{
    A = lah_matLoad("./test_matrices/TestMatrix_matAdd1.dat");
    B = lah_matLoad("./test_matrices/TestMatrix_matAdd2.dat");
    C_exp = lah_matLoad("./test_matrices/SolMatrix_matAdd.dat");

    C_res = lah_matAlloc(A->nC, A->nR, 1);
}

void test_teardown(void) 
{
    lah_matFree(A);
    lah_matFree(B);
    lah_matFree(C_res);
    lah_matFree(C_exp);
}

MU_TEST(add_a_b)
{
    mu_assert(!lah_matAdd(C_res, 1 , A, 1, B), "Error when calling lah_matAdd()");
    mu_assert_mat_eq(C_exp, C_res, 4);
}

MU_TEST(check_2A) 
{
    mu_assert(!lah_matAdd(B, 1 , A, 1, A), "Error when calling lah_matAdd()");
    mu_assert(!lah_matAdd(C_res, 2 , A, 0, A), "Error when calling lah_matAdd()");
    
    mu_assert_mat_eq(B, C_res, 4);
}

MU_TEST_SUITE(test_suite) {
	MU_SUITE_CONFIGURE(&test_setup, &test_teardown);
	MU_RUN_TEST(add_a_b);
	MU_RUN_TEST(check_2A);
}

int main() 
{
    MU_RUN_SUITE(test_suite);
	
    MU_REPORT();

	return 0;
}
