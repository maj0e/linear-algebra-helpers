/*****************************************************************************/
/* Solvers for the Least squares problem min_x || B - A*x ||                 */
/*****************************************************************************/
#include <lah.h>

/* lah_solveLS
 * ------------
 * Solves the least squares problem by applying the Lapack routine gels
 * which uses QR or LQ factorization.
 * Params:  - A, B, X: Matrices of the LS problem min_x || B - A*x ||
*/
#ifdef HAVE_LAPACK /* Use a LAPACK package to do the heavy work */
#include <lapacke.h>
lah_Return lah_solveLS(koop_mat *A, *B, *X)
{
    lah result = 0;
    gels('N', A->nR, A->nC, NRHS, A->data, LDA, B, LDB, WORK,
                        LWORK, result);
    
    /* Alternative: used for ill conditioned problems
     * ca. 5 times slower. Could check residuum output to determine which precision to use*/
    /* gelsd () */
    return result;
}
#endif
