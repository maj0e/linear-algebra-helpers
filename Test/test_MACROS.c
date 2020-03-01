#include "lah.h"
#include "minunit.h"

lah_mat *C_exp;
lah_mat *C_res;
lah_mat *A;
lah_mat *B;

/* Number of significant digits to compare floats*/
#define PRECISION 4

MU_TEST(Squared)
{
    A = lah_matAlloc(2, 2, 1);
    B = lah_matAlloc(1, 2, 1);
    mu_check(A != NULL);
    mu_check(B != NULL);
    
    mu_check(LAH_ISTYPE(A, lahTypeSquare));
    mu_check(!LAH_ISTYPE(B, lahTypeSquare));
    
    mu_check(NULL == lah_matFree(A));
    mu_check(NULL == lah_matFree(B));
}

MU_TEST(NoData) 
{
    A = lah_matAlloc(2, 2, 0);
    B = lah_matAlloc(1, 2, 1);
    mu_check(A != NULL);
    mu_check(B != NULL);
    
    mu_check(LAH_ISTYPE(A, lahTypeNoData));
    mu_check(!LAH_ISTYPE(B, lahTypeNoData));
    
    mu_check(NULL == lah_matFree(A));
    mu_check(NULL == lah_matFree(B));
}

MU_TEST(MultipleTypes) 
{
    A = lah_matAlloc(1, 1, 0);
    mu_check(A != NULL);
    
    mu_check(LAH_ISTYPE(A, lahTypeNoData));
    mu_check(LAH_ISTYPE(A, lahTypeSquare));
    mu_check(LAH_ISTYPE(A, lahTypeVector));
    
    mu_check(NULL == lah_matFree(A));
}

MU_TEST(Set_Unset_Types)
{
    A = lah_matAlloc(2, 2, 0);
    mu_check(A != NULL);
    
    LAH_SETTYPE(A, lahTypePositive);
    LAH_SETTYPE(A, lahTypeSymmetric);
    LAH_SETTYPE(A, lahTypeDiagonal);
    LAH_SETTYPE(A, lahTypeOrthogonal);
    
    mu_check(LAH_ISTYPE(A, lahTypePositive));
    mu_check(LAH_ISTYPE(A, lahTypeSymmetric));
    mu_check(LAH_ISTYPE(A, lahTypeDiagonal));
    mu_check(LAH_ISTYPE(A, lahTypeOrthogonal));
    
    LAH_UNSETTYPE(A, lahTypePositive);
    LAH_UNSETTYPE(A, lahTypeSymmetric);
    
    
    mu_check(!LAH_ISTYPE(A, lahTypePositive));
    mu_check(!LAH_ISTYPE(A, lahTypeSymmetric));
    mu_check(LAH_ISTYPE(A, lahTypeDiagonal));
    mu_check(LAH_ISTYPE(A, lahTypeOrthogonal));
    
    LAH_UNSETTYPE(A, lahTypeDiagonal);
    LAH_UNSETTYPE(A, lahTypeOrthogonal);
    
    mu_check(! LAH_ISTYPE(A, lahTypePositive));
    mu_check(! LAH_ISTYPE(A, lahTypeSymmetric));
    mu_check(! LAH_ISTYPE(A, lahTypeDiagonal));
    mu_check(! LAH_ISTYPE(A, lahTypeOrthogonal));
    
    mu_check(NULL == lah_matFree(A));
}

MU_TEST(Toggle_Types)
{
    A = lah_matAlloc(2, 2, 0);
    mu_check(A != NULL);
    
    LAH_SETTYPE(A, lahTypeSymmetric);
    LAH_TOGGLETYPE(A, lahTypeSymmetric);
    mu_check(!LAH_ISTYPE(A, lahTypeSymmetric));
    LAH_TOGGLETYPE(A, lahTypeSymmetric);
    mu_check(LAH_ISTYPE(A, lahTypeSymmetric));
    
    mu_check(NULL == lah_matFree(A));
}
    
MU_TEST(Mat_Entries)
{
    A = lah_matAlloc(2, 2, 1);
    A->data[0] = 1;
    A->data[A->incRow] = 2;
    A->data[A->incCol] = 3;
    A->data[A->incCol + A->incRow] = 4;
    
    mu_check(1 == LAH_ENTRY(A, 0, 0));
    mu_check(2 == LAH_ENTRY(A, 1, 0));
    mu_check(3 == LAH_ENTRY(A, 0, 1));
    mu_check(4 == LAH_ENTRY(A, 1, 1));
    
    mu_check(NULL == lah_matFree(A));
}


MU_TEST_SUITE(test_suite) {
    MU_RUN_TEST(Squared);
    MU_RUN_TEST(NoData);
    MU_RUN_TEST(MultipleTypes);
    MU_RUN_TEST(Set_Unset_Types);
    MU_RUN_TEST(Toggle_Types);
	MU_RUN_TEST(Mat_Entries);
}

int main() 
{
    MU_RUN_SUITE(test_suite);
	
    MU_REPORT();

	return 0;
}

