#include "lah.h"

#include <math.h>
#include <stdlib.h>

/* Needed for some examples*/
lah_value lah_gaussNoise()
{
    lah_value val1, val2, w;
    do
    {
        val1 = (lah_value)(2.0 * rand())  / RAND_MAX - 1.0;
        val2 = (lah_value)(2.0 * rand())  / RAND_MAX - 1.0;
        w = val1 * val1 + val2 * val2;
    } while (w >= 1.0);

    return val1 * sqrt((-2.0 * log(w)) / w);
}

/* Maybe there is also a more sophisticated method
 * like ther Marsaglia algorithm for Gaussian Noise */
lah_value lah_lorentzNoise()
{
     return tan(LAH_PI * (lah_value)(rand()) / RAND_MAX - 0.5);
}
