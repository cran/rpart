/*  SCCS @(#)formatg.c	1.1 03/31/01 */
/*
** Write a bunch of numbers in a desired C format, (almost always %g)
*/

#include <stdio.h>
#include "rpartS.h"

void formatg(Sint *n, double *x, char **format, char **out) {
    int i;

    for (i=0; i< *n; i++) {
	sprintf(out[i], format[i], x[i]);
	}
    }
