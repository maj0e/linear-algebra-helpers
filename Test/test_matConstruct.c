/*******************************************************************/
/** Test of the Cholesky factorization                            **/
/** Should be tested with all combinations of HAVE_LAPACK,        **/
/** Row-/ColMajor and USE_SINGLE                                  **/
/*******************************************************************/

#include "lah.h"
#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PRECISION 5


#define MAX_VALUE 100.0

#define NUMBER_TESTS 6
lah_index test_size_array[NUMBER_TESTS] = {1, 2, 5, 10, 20, 30};
lah_index test_size = 0;
lah_index nParameters = 0;

lah_value *parameter_vec;
lah_mat *Parameter_Mat;
lah_mat *RandMat;

void test_setup(void) 
{
    lah_index i;
    nParameters = (test_size * (test_size + 3)) / 2;
    parameter_vec = malloc(nParameters * sizeof(lah_value));
    srand((unsigned)time(NULL));
    for (i = 0; i < nParameters ; i++)
    {
        parameter_vec[i] = ((lah_value)rand() / RAND_MAX - 0.5)*2.0*MAX_VALUE;
    }
    Parameter_Mat = lah_matAlloc(test_size, test_size, 1);
    for (i = 0; i < test_size * test_size ; i++)
    {
        Parameter_Mat->data[i] = ((lah_value)rand() / RAND_MAX - 0.5)*2.0*MAX_VALUE;
    }
}

void test_teardown(void) 
{
    free(parameter_vec);
    lah_matFree(RandMat);
    lah_matFree(Parameter_Mat);
}

MU_TEST(test_construct_random_positive)
{
    RandMat = lah_matConstructPo(test_size, parameter_vec);
    mu_check(RandMat != NULL);
    mu_check(lah_isPositive(RandMat));
}

MU_TEST(test_construct_random_positive_simple)
{
    RandMat = lah_matConstructPo_simple(Parameter_Mat);
    mu_check(RandMat != NULL);
    mu_check(lah_isPositive(RandMat));
}

MU_TEST_SUITE(test_constructor_suite)
{   
	lah_index i;
	MU_SUITE_CONFIGURE(&test_setup, &test_teardown); 
	for(i = 0; i < NUMBER_TESTS; i++)
	{
        test_size = test_size_array[i];
		MU_RUN_TEST(test_construct_random_positive);
        MU_RUN_TEST(test_construct_random_positive_simple);
	}
}
/*******************************************************************/

int main()
{
    MU_RUN_SUITE(test_constructor_suite);
    MU_REPORT();
    
    return 0;
}
