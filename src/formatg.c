/*  SCCS @(#)formatg.c	1.2 06/06/01 */
/*
** Write a bunch of numbers in a desired C format, (almost always %g)
*/
#include <stdio.h>
#include "rpartS.h"
#ifdef WIN32
#include <ctype.h>
#endif
void formatg( Sint *n, double *x, char **format, char **out) 
{
    int i;
#ifdef WIN32
    int len;
    char *p;
#endif

    for (i=0; i< *n; i++) {
	sprintf(out[i], format[i], x[i]);
#ifdef WIN32
	/* change e+/-00n to e+/-0n etc */
	p = out[i];
	len = strlen(p);
	if (tolower(p[len-5]) == 'e' &&
	    (p[len-4] == '+' || p[len-4] == '-') &&
	    p[len-3] == '0' &&
	    isdigit(p[len-2]) && isdigit(p[len-1])) {
	    p[len-3] = p[len-2];
	    p[len-2] = p[len-1];
	    p[len-1] = '\0';
	}
#endif
    }
}
