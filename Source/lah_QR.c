#include "lah.h"

#ifdef HAVE_LAPACK /* Use a LAPACK package to do the heavy work */
#include <lapacke.h>
lah_Return lah_QR(lah_mat *A, lah_value *tau)
{
    if (A == NULL || tau == NULL) 
        return lahReturnParameterError;
    
    /* Needs leading dimensions */ 
    /*lda = (matrix_layout == LAPACK_COL_MAJOR) ? A->nR : A->nC;
    */    
    return (!GEQRF(LAH_LAPACK_LAYOUT, A->nR, A->nC, A->data , LAH_LEADING_DIM(A), tau)) 
    ? lahReturnOk : lahReturnExternError;
}

#else /* fallback to primitive linear algebra subroutines */
/* Copied from rosetta code, use as  a blueprint*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
typedef struct {
	int m, n;
	double ** v;
} mat_t, *mat;
 
mat matrix_new(int m, int n)
{
	mat x = malloc(sizeof(mat_t));
	x->v = malloc(sizeof(double*) * m);
	x->v[0] = calloc(sizeof(double), m * n);
	for (int i = 0; i < m; i++)
		x->v[i] = x->v[0] + n * i;
	x->m = m;
	x->n = n;
	return x;
}
 
void matrix_delete(mat m)
{
	free(m->v[0]);
	free(m->v);
	free(m);
}
 
void matrix_transpose(mat m)
{
	for (int i = 0; i < m->m; i++) {
		for (int j = 0; j < i; j++) {
			double t = m->v[i][j];
			m->v[i][j] = m->v[j][i];
			m->v[j][i] = t;
		}
	}
}
 
mat matrix_copy(int n, double a[][n], int m)
{
	mat x = matrix_new(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			x->v[i][j] = a[i][j];
	return x;
}
 
mat matrix_mul(mat x, mat y)
{
	if (x->n != y->m) return 0;
	mat r = matrix_new(x->m, y->n);
	for (int i = 0; i < x->m; i++)
		for (int j = 0; j < y->n; j++)
			for (int k = 0; k < x->n; k++)
				r->v[i][j] += x->v[i][k] * y->v[k][j];
	return r;
}
 
mat matrix_minor(mat x, int d)
{
	mat m = matrix_new(x->m, x->n);
	for (int i = 0; i < d; i++)
		m->v[i][i] = 1;
	for (int i = d; i < x->m; i++)
		for (int j = d; j < x->n; j++)
			m->v[i][j] = x->v[i][j];
	return m;
}
 
/* c = a + b * s */
double *vmadd(double a[], double b[], double s, double c[], int n)
{
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return c;
}
 
/* m = I - v v^T */
mat vmul(double v[], int n)
{
	mat x = matrix_new(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			x->v[i][j] = -2 *  v[i] * v[j];
	for (int i = 0; i < n; i++)
		x->v[i][i] += 1;
 
	return x;
}
 
/* ||x|| */
double vnorm(double x[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
}
 
/* y = x / d */
double* vdiv(double x[], double d, double y[], int n)
{
	for (int i = 0; i < n; i++) y[i] = x[i] / d;
	return y;
}
 
/* take c-th column of m, put in v */
double* mcol(mat m, double *v, int c)
{
	for (int i = 0; i < m->m; i++)
		v[i] = m->v[i][c];
	return v;
}
 
void matrix_show(mat m)
{
	for(int i = 0; i < m->m; i++) {
		for (int j = 0; j < m->n; j++) {
			printf(" %8.3f", m->v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
 
double in[][3] = {
	{ 12, -51,   4},
	{  6, 167, -68},
	{ -4,  24, -41},
	{ -1, 1, 0},
	{ 2, 0, 3},
};
 
int main()
{
	mat R, Q;
	mat x = matrix_copy(3, in, 5);
	householder(x, &R, &Q);
 
	puts("Q"); matrix_show(Q);
	puts("R"); matrix_show(R);
 
	// to show their product is the input matrix
	mat m = matrix_mul(Q, R);
	puts("Q * R"); matrix_show(m);
 
	matrix_delete(x);
	matrix_delete(R);
	matrix_delete(Q);
	matrix_delete(m);
	return 0;
}



#endif /* endif primitive implementation */
