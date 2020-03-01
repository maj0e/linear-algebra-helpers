#include "minunit.h"
#include "lah.h"
#include <stdio.h>

#define PRECISION -1

lah_mat *P;
lah_mat *Psqrt;
lah_mat *Work;

void test_setup(void) 
{
    P = lah_matLoad("test_matrices/TestMatrix_matSqrt.dat");
    Psqrt = lah_matLoad("test_matrices/TestMatrix_matSqrt.dat");
    Work = lah_matAlloc(P->nC, P->nR, 1);
}

void test_teardown(void) 
{
    lah_matFree(P);
    lah_matFree(Psqrt);
    lah_matFree(Work);
}

MU_TEST(Psqrt_squared_eq_P)
{
    mu_check(!lah_matSqrt(P, Psqrt, Work));
    lah_matMul(lahNorm, lahNorm, 0.0, 1.0, Work, Psqrt, Psqrt);
    mu_assert_mat_eq(P, Work, PRECISION);
}

MU_TEST_SUITE(matSqrt_test_suite)
{
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	MU_RUN_TEST(Psqrt_squared_eq_P);
}

int main ()
{
    MU_RUN_SUITE(matSqrt_test_suite);
    MU_REPORT();
    
    return 0;
}
