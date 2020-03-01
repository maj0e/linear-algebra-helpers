#include "lah.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>

/* Construct a diagonal matrix from n parameters */
lah_mat *lah_matConstructDiag(lah_index nC, const lah_value *parameter_vec)
{
	lah_index i, j;
    lah_index counter = 0;
    lah_index isPositive = 1;
	lah_mat *Res = lah_matAlloc(nC, nC, 1);
    if (parameter_vec == NULL)
    {
        for(i = 0; i < Res->nR; i++)
        {
            for(j = 0; j < Res->nC; j++)
            {
                LAH_ENTRY(Res, i, j) = 0;
            }
            LAH_ENTRY(Res, i, i) = 1;
        }
    }
    else
    {
        for(i = 0; i < Res->nC; i++)
        {
            LAH_ENTRY(Res, i, i) = parameter_vec[counter];
            if (parameter_vec[counter] < 0)
                isPositive = 0;
            counter++;
        }	
        
    }
	LAH_SETTYPE(Res, lahTypeSquare);
	LAH_SETTYPE(Res, lahTypeDiagonal);
    if(isPositive == 1)
        LAH_SETTYPE(Res, lahTypeDiagonal);
	return Res;
}

/*
lah_mat *lah_matGeneralToSymmetric(lah_mat *A)
{
	lah_mat *Res = lah_matAlloc(A->nC, A->nR, 1);
	lah_mat *At = lah_matGetTransView(A);
	lah_matAdd(A, 0.5, A, 0.5, At);
	return Res;
}
*/

/* Construct a symmetric Matrix from n(n+1)/2 parameters */
lah_mat *lah_matConstructSy(lah_index nC, const lah_value *parameter_vec)
{
	lah_index i, j;
    lah_index counter = 0;
	lah_mat *Res = lah_matAlloc(nC, nC, 1);
	for(i = 0; i < Res->nC; i++)
	{
		LAH_ENTRY(Res, i, i) = parameter_vec[counter];
		counter++;
		for(j = i + 1; j < Res->nR; j++)
		{
			LAH_ENTRY(Res, j, i) = parameter_vec[counter];
			LAH_ENTRY(Res, i, j) = parameter_vec[counter];
			counter++;
		}
	}	
	LAH_SETTYPE(Res, lahTypeSquare);
	LAH_SETTYPE(Res, lahTypePositive);
	return Res;
}

/* Construct a orthogonal Matrix from n(n+1)/2 parameters */
lah_mat *lah_matConstructOrtho(lah_index nC, lah_value *parameter_vec)
{
    /* Subgroup Algorithm */
    lah_index j, xlength;
    int *signs;
    int k;
    lah_value norm;
    lah_value beta;
    lah_mat *x;
    lah_mat *y;
    lah_mat *A;
    lah_mat *A_sub;
    
    signs = malloc(nC * sizeof(int));
    A = lah_matConstructDiag(nC, NULL); /* Identity Matrix*/
    y = lah_matAlloc(nC, 1, 1);
    x = lah_matAlloc(1, nC, 0);
    A_sub = lah_matAlloc(nC, nC, 0);
    x->data = parameter_vec;
    /* set last row to random sign*/
    if (x->data[0] < 0)
        signs[nC-1] = -1;
    else
        signs[nC-1] = 1;
    x->data++;
    for (k = nC-2; k >= 0; k--)//(k = n-1 to 1 by -1; 
    {    
        /* generate random Householder transformation */
        /*Create Random vector with 1 column and n-k+1 rows*/
        xlength = nC-k;
        x->nR = xlength;
        
        norm = 0;
        for (j = 0; j < xlength; j++)
        {
            norm += x->data[j] * x->data[j];
        }
        norm = sqrt(norm);
        
        if (x->data[0] < 0)
        {
            signs[k] = 1;
            norm = (-1) * norm;
        }
        else
            signs[k] = -1;
        
        x->data[0] += norm;                  /* Add norm to first entry*/
        beta = 1.0 / (norm * x->data[0]);    /* beta = s(x[1] + s) */
        A_sub->data = A->data + k * (A->incRow + A->incCol);
        A_sub->nC = xlength;
        A_sub->nR = xlength;
        /* apply the transformation to A */
        y->nC = xlength;
        lah_matMul(lahTrans, lahNorm, 0.0, 1.0, y, x, A_sub); /* y = x^t*A[k:n, ] */
        lah_matMul(lahNorm, lahNorm, 1.0, beta, A_sub, x, y); /* A[k:n, ] = A[k:n, ] - beta * x * y */
        x->data += xlength; /* Use the next parameters */
    }
    
    /* change signs of k_th row when signs[k]=-1 */
    for (k = nC-1; k >= 0; k--)
    {   
        if (signs[k] < 0)
        {
            for(j = 0; j < nC; j++)
            {
                LAH_ENTRY(A, k, j) = (-1) * LAH_ENTRY(A, k, j);
            }
        }
        
    }
    free(signs);
    free(x);
    free(A_sub);
    lah_matFree(y);
    return A;
}
lah_mat *lah_matConstructPo_simple(lah_mat *A)
{
    lah_mat *P;
    if (!LAH_ISTYPE(A, lahTypeSquare))
        return NULL;
    
    /* A = A * A^T*/
    P = lah_matAlloc(A->nC, A->nR, 1);
    lah_matMul(lahNorm, lahTrans, 0.0, 1.0, P, A, A);
    return P;
}

lah_mat *lah_matConstructPo(lah_index nC, lah_value *parameter_vec)
{
	/* Construct a positive definite matrix by P = Q D Q^t */
 	/* Q is an orthogonal matrix, which can be constructed by Householder reflections */
 	/* D is a diagonal matrix holding the EigenValues */
 	lah_index i, j;
	lah_mat *Res = lah_matAlloc(nC, nC, 1);
    lah_mat *Work = lah_matAlloc(nC, nC, 1);	
	lah_mat *Q = lah_matConstructOrtho(nC, parameter_vec + nC);
	/* Calculate Res = D Q^T*/
	for(i = 0; i < Q->nC; i++)
	{
		for(j = 0; j < Q->nR; j++)
		{
			LAH_ENTRY(Work, i, j) = fabs(parameter_vec[i] +0.1) * LAH_ENTRY(Q, j, i);	
		}
	}
    /* Caclulate Res = Q * Work = Q * D * Q^t */
	lah_matMul(lahNorm, lahNorm, 0.0, 1.0, Res, Q, Work);
	lah_matFree(Q);
    lah_matFree(Work);
	LAH_SETTYPE(Res, lahTypeSquare);
	LAH_SETTYPE(Res, lahTypeSymmetric);
	LAH_SETTYPE(Res, lahTypePositive);
	return Res; 
}
