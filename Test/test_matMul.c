#include "lah.h"
#include "minunit.h"

lah_mat *C_exp;
lah_mat *C_res;
lah_mat *A;
lah_mat *B;
lah_mat *D;

/* Number of significant digits to compare floats*/
#define PRECISION 4

void test_setup(void) 
{
    A = lah_matLoad("./test_matrices/TestMatrix_matMul1.dat");
    B = lah_matLoad("./test_matrices/TestMatrix_matMul2.dat");
    C_exp = lah_matLoad("./test_matrices/SolMatrix_matMul.dat");
    if (A ==NULL || B == NULL || C_exp == NULL)
        mu_fail("Matrices not loaded!!!\n");
    C_res = lah_matAlloc(A->nC, A->nR, 1);
}

void test_teardown(void) 
{
    lah_matFree(A);
    lah_matFree(B);
    lah_matFree(C_res);
    lah_matFree(C_exp);
}

MU_TEST(C_eq_A_x_B)
{
    mu_assert(!lah_matMul(lahNorm, lahNorm, 0.0, 1.0, C_res, A, B), "Error when calling lah_matMul()");
    mu_assert_mat_eq(C_exp, C_res, PRECISION);
}

MU_TEST(Ct_eg_Bt_x_At) 
{
    lah_mat *Ct_exp = lah_matTrans(C_exp);
    mu_assert(!lah_matMul(lahTrans, lahTrans, 0.0, 1.0, C_res, B, A), "Error when calling lah_matMul()");
    
    mu_assert_mat_eq(Ct_exp, C_res, PRECISION);
    lah_matFree(Ct_exp);
}

MU_TEST_SUITE(test_suite) {
	MU_SUITE_CONFIGURE(&test_setup, &test_teardown);
	MU_RUN_TEST(C_eq_A_x_B);
	MU_RUN_TEST(Ct_eg_Bt_x_At);
}

int main() 
{
    MU_RUN_SUITE(test_suite);
    MU_REPORT();

	return 0;
}
