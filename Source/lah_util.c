#include "lah.h"

#include <stdlib.h>
#include <stdio.h>

/* allocation for the Matric struct lah_mat */
lah_mat *lah_matAlloc(lah_index nC, lah_index nR,
                      lah_index AllocData)
{
    /* Allocate the struct */
    lah_mat *mat = NULL;
    size_t size = 0;
    mat = (lah_mat *) calloc(1, sizeof(lah_mat));
    if (mat == NULL) return (NULL) ;                /* out of memory */
    mat->nC = nC ;                       /* define dimensions */
    mat->nR = nR ;
    
    mat->matType = 0; /* Should already be zero */
    
    /* Vector ? */
    if (nC == 1)
        LAH_SETTYPE(mat, lahTypeVector);
    
    /* Square-Matrix? */
    if (nC == nR)
        LAH_SETTYPE(mat, lahTypeSquare);

    if (LAH_LAYOUT == lahRowMajor)
    {
        mat->incRow = nC; 
        mat->incCol = 1;
        LAH_UNSETTYPE(mat, lahTypeColMajor);
    }
    else
    {
        mat->incRow = 1;
        mat->incCol = nR; 
        LAH_SETTYPE(mat, lahTypeColMajor);
    }
    /* if AllocData == 1 allocate place for the data, 
     * otherwise set null, so that the pointer can be set later */
    if (AllocData != 0)
    {
        size = nC * nR;
        mat->data = (lah_value*) calloc(size, sizeof(lah_value));
    }
    else
    {
        mat->data = NULL;
    }
    
    /* Owns data ? */
    if (AllocData == 0)
        LAH_SETTYPE(mat, lahTypeNoData);
    
    if( (AllocData == 1) && (mat->data == NULL))
        return lah_matFree(mat);
    
    return mat;
}

/* free Matrix struct lah_mat */
lah_mat *lah_matFree(lah_mat *mat)
{
    if (!mat) return NULL;      /* do nothing if A already NULL */
    if(!LAH_ISTYPE(mat, lahTypeNoData))
        free(mat->data);
    free(mat);
    return NULL; 
}

/* load a matrix from a file */ 
lah_mat *lah_matLoad(const char *fname)
{
    FILE *file = NULL;
    /* use auxilary variable, because type of 
     * lah_value ist not exactly known */
    lah_value temp = 0.0;
    lah_index nC = 0;
    lah_index nR = 0;
    lah_index i, j;
    lah_index counter = 0;
    lah_mat *mat = NULL;
    
    file = fopen(fname, "r");
    if (file == NULL) return NULL;
    if (1 != fscanf(file, "%lu", &nC))
        return NULL;
    /*nC = (lah_index) temp;*/
    
    if (1 != fscanf(file, "%lu", &nR))
        return NULL;
    /*nR = (lah_index) temp; */

    mat = lah_matAlloc(nC, nR, 1);
    if (!mat) return NULL;
    
    /* Does not enforce order of data */
    for (j = 0; j < nR; j++)
    {
        for (i = 0; i < nC; i++)
        {
            counter += fscanf(file, "%lg", &temp);
            LAH_ENTRY(mat, j, i) = (lah_value) temp;
        }
    }
    fclose(file);
    
    if (counter != nC * nR)
        return lah_matFree(mat);
    else
        return mat;
}

/* print a lah_mat; use %g for integers to avoid differences 
 * with lah_index */
lah_Return lah_matPrint (lah_mat const *A, lah_index brief)
{
    lah_index i, j;
    /*lah_index incRow = A->incRow;
    lah_index incCol = A->incCol;
*/
    if (!A) 
    { 
        printf ("(null)\n") ; 
        return (lahReturnParameterError); 
    }
    
    printf ("%g-by-%g, incRow: %g incCol: %g\n", (double) A->nR,
            (double) A->nC, (double) A->incRow, (double) A->incCol);

    for (i = 0 ; i < A->nR ; i++)
    {
       for (j = 0 ; j < A->nC ; j++)
       {
           printf("  %g,", (lah_value) *(A->data + A->incRow * i + A->incCol * j)) ;
           if (brief && j > 10) 
           {
               printf ("  ...\n");
               break; 
           }
       }
       printf("\n");
       if (brief && i > 10) 
       { 
           printf ("...\n...\n"); 
           return lahReturnOk;
       }
    }
    printf("\n");
    return lahReturnOk ;
}
