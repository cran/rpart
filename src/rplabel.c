/* SCCS @(#)rplabel.c	1.6 02/08/98 */
/*
** Called by S to make nice "labels" for the split points
**
**  csplit is a matrix with maxcat rows and ncat columns
**  splits has (nsplit-1)/2 columns, 2 rows
**
**  cutleft and cutright are the output
*/
#include "rpart.h"
#include "rpartproto.h"
#include <stdlib.h>
#include <stdio.h>
static char *Rp_strdup();

void rplabel(int *nsplit,   int *index,   double *splits, 
             int *ncat,     int *csplit,  char   **cutleft, char **cutright)
    {
    int i,j, k;
    int  ii, jj;
    int  nn, nc, ccol;
    char buf[1000];

    nn = (*nsplit-1) /2;    /*leaves = splits +1;  *nsplit = leaves + splits*/
    nc = *ncat;
    j=0;
    for (i=0; i< *nsplit; i++) {
	if (index[i]==0) continue;

	k = splits[j];
	if (k >=2)  {  /* categorical */
	    ccol = splits[j + nn] -1;
	    buf[0] = ':';
	    ii =1;
	    for (jj=0; jj< k; jj++)
		if      (csplit[jj*nc +ccol] == LEFT)  buf[ii++] = 'a' + jj;
	    buf[ii] = '\0';
	    cutleft[i] = Rp_strdup(buf);

	    buf[0] = ':';
	    ii =1;
	    for (jj=0; jj< k; jj++)
		if      (csplit[jj*nc + ccol] == RIGHT)  buf[ii++] = 'a' + jj;
	    buf[ii] = '\0';
	    cutright[i] = Rp_strdup(buf);
	    }

	else {  /* continuous */
	    if (k == LEFT) {
		sprintf(buf, "<%.6g", splits[nn +j]);
		cutleft[i] = Rp_strdup(buf);
		sprintf(buf, ">%.6g", splits[nn +j]);
		cutright[i] = Rp_strdup(buf);
		}
	    else {
		sprintf(buf, ">%.6g", splits[nn +j]);
		cutleft[i] = Rp_strdup(buf);
		sprintf(buf, "<%.6g", splits[nn +j]);
		cutright[i] = Rp_strdup(buf);
		}
	    }
	j++;
	}
    }

#include "rpartS.h"
static char *Rp_strdup(s)
char *s;
{
	return(strcpy(CALLOC(strlen(s)+1, 1), s));
}
