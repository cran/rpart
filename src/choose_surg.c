/* SCCS @(#)choose_surg.c	1.3 02/08/98 */
/*
** A particular split routine, optimized for the surrogate variable
**  search.  The "goodness" of a split is the total number of concordant
**  observations between the surrogate and the primary split.
**  Note that the CART folks use the %concordance, which factors missing
**  values into the equations somewhat differently.
**
**  y is coded as  +1=left, -1=right, 0=missing 
*/
#include "rpart.h"
#include "rpartproto.h"

void choose_surg(int nodenum,    int *y,         double *x,     long *order, 
		 int ncat,       int *agreement, double *split, int *csplit)
    {
    int i,j;
    int agree;
    int lcount, rcount;
    int ll, lr, rr, rl;
    int defdir;
    double lastx;
    int  *which, *left, *right;

    which = rp.which;
    left = rp.left;
    right =rp.right;

    if (ncat==0) {  /* continuous case */
	/*
	** ll = y's that go left that are also sent left by my split
	** lr = y's that go left that I send right
	** rl= y's that go right that I send to the left
	** rr= y's that go right that I send to the right
	**
	** The agreement is max(ll+rr, lr+rl)
	**
	** I enforce that at least 2 obs must go each way, to avoid having an
	**  uncorrelated surrogate beat the "null" surrogate too easily
	*/
	ll = rl =0;
	for (i=rp.n-1; i>=0; i--) { /*start with me sending all to the left */
	    j = order[i];
	    if (j>=0 &&  which[j]==nodenum) {
		lastx = x[i];        /*this is why I ran the loop backwards*/
		switch( y[j]) {
		    case LEFT : ll++;
				break;
		    case RIGHT: rl++;
				break;
		    default:;
		    }
		}
	    }
	lr = rr =0;
	if (ll>rl) agree = ll;
	else       agree = rl;

	for (i=0; (ll+rl)>=2; i++) {
	    j = order[i];
	    if (j >=0 &&  which[j]==nodenum) {
		if ((lr+rr)>=2  &&  x[i] != lastx) {
		    if ((ll+rr) > agree) {
			agree = ll + rr;
			csplit[0] = RIGHT;       /* < goes to the right */
			*split = (x[i] + lastx)/2;
			}
		    else if ((lr+rl) > agree) {
			agree = lr + rl;
			csplit[0] = LEFT;
			*split = (x[i] + lastx)/2;
			}
		    }

		switch (y[j]) {
		    case LEFT : ll--; lr++;
				break;
		    case RIGHT: rl--; rr++;
				break;
		    default: ;         /*ignore missing y's */
		    }
		lastx = x[i];
		}
	    }
	}

    else {     /* categorical predictor */
	for (i=0; i<ncat; i++) {
	    left[i] =0;
	    right[i]=0;
	    }
	for (i=0; i<rp.n; i++) {
	    if (which[i] == nodenum &&  order[i]>=0) {
		j = x[i] -1;
		switch( y[i]) {
		    case LEFT : left[j]++;
				break;
		    case RIGHT: right[j]++;
				break;
		    default:;
		    }
		}
	    }

	/*
	**  Tied categories get sent to the majority, if possible
	*/
	lcount=0; rcount=0;
	for (i=0; i<ncat; i++) {
	    lcount += left[i];
	    rcount += right[i];
	    }
	if (lcount > rcount) defdir = LEFT;
	else                 defdir = RIGHT;

	agree=0;  /*agreement */
	for (i=0; i<ncat; i++) {
	    if (left[i]==0 && right[i]==0) csplit[i]=0;
	    else {
		if (left[i]< right[i] || (left[i]==right[i] && defdir==RIGHT)) {
		    agree+= right[i];
		    csplit[i] = RIGHT;
		    }
		else {
		    agree += left[i];
		    csplit[i] = LEFT;
		    }
		}
	    }
	}
    *agreement = agree;
    }
