#include "lah.h"

/* This will return a lah_mat struct At with inverted order 
 * (RowMajor <-> ColMajor) */
lah_mat *lah_matTrans(const lah_mat *A)
{
    lah_mat *At = lah_matAlloc(A->nR, A->nC, 0);

    /*
    if (LAH_LAYOUT == lahRowMajor)
    {
        At->incRow = 1;
        At->incCol = A->nR; 
    }
    else
    {
        At->incRow = A->nC; 
        At->incCol = 1;
    }*/
    At->incRow = A->incCol;
    At->incCol = A->incRow;
    
    At->data = A->data;
    return At;
}

lah_mat *lah_matTrans_copy(lah_mat *A)
{
    lah_index i, j;
    lah_mat *At = lah_matAlloc(A->nC, A->nR, 1);
    
    for (i = 0; i < A->nR; i++)
    {
        for (j = 0; j < A->nC; j++)
        {
            LAH_ENTRY(At, j, i) = LAH_ENTRY(At, i, j);
        }
    }
    
    return At;
    
}
