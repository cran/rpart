/* SCCS @(#)poissonpred.c	1.3  02/08/98 */
/*
** return the additive component of the deviance for a particular point
*/
#include <math.h>

double poissonpred(double *y, double *lambda)
    {
    double temp, dev;

    temp = y[1];
    dev  = temp - *lambda*y[0];
    if (temp>0)
	dev += temp * log(*lambda * y[0] / temp);

    return(-2*dev);
    }
