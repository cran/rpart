/* SCCS @(#)surrogate.c	1.5 02/08/98 
/*
** Calculate the surrogate splits for a node and its primary
**    (This routine is an awful lot like bsplit)
**
** Input :      node
**              start and stop indices for the arrays (which obs apply)
**
** Output:      Fills in the node's
**                      surrogate splits
**                      lastsurrogate value
**
** Uses:        The global vector tempvec (integer) as a temporary, assumed
**                to be of length n.
*/
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

void surrogate(struct node *me, int nodenum)
    {
    int i, j, k;
    int var;   /* the primary split variable */
    int agreement;         /*number of agreements in the split */
    double split;
    double improve;
    int lcount, rcount;    /* number sent left and right by primary */
    int count;
    unsigned int extra;
    struct split *ss;
    int  *index;
    int  *which,
	 *tempy;
    long **sorts;
    double **xdata;
    int ncat;

    which = rp.which;
    tempy = rp.tempvec;
    sorts = rp.sorts;
    xdata = rp.xdata;
    /*
    ** First construct, in tempy, the "y" variable for this calculation.
    ** It will be LEFT:goes left, 0:missing, RIGHT:goes right.
    **  Count up the number of obs the primary sends to the left, as my
    **  last surrogate (or to the right, if larger).
    */
    var = (me->primary)->var_num;
    if (rp.numcat[var]==0)  {  /* continuous variable */
	split = (me->primary)->spoint;
	extra = (me->primary)->csplit[0];
	for (i=0; i<rp.n; i++) {
	    j = sorts[var][i];
	    if (j<0) tempy[-(j+1)]=0;
	    else if (which[j] == nodenum) {
		if (xdata[var][i] < split)
			 tempy[j] = extra;
		else
			 tempy[j] =  -extra;
		}
	    }
	}

    else {  /* categorical variable */
	index = (me->primary)->csplit;
	for (i=0; i<rp.n; i++) {
	    if (which[i] != nodenum) continue;
	    if (sorts[var][i]<0) tempy[i] =0;
	    else        tempy[i] = index[(int)xdata[var][i] -1];
	    }
	}

    lcount=0; rcount=0;
    for (i=0; i<rp.n; i++) {
	if (which[i] != nodenum) continue;
	switch(tempy[i]) {
	    case LEFT : lcount++;  break;
	    case RIGHT: rcount++;  break;
	    default: break;
	    }
	}

    me->lastsurrogate = lcount;
    if (lcount < rcount) count = rcount;
      else               count = lcount;
    /* surrogates that don't get a least "count" correct aren't kept */

    /*
    ** Now walk through the variables
    */
    me->surrogate =0;
    for (i=0; i<rp.nvar; i++) {
	if (var == i) continue;
	ncat = rp.numcat[i];

	choose_surg(nodenum, tempy, xdata[i], sorts[i], ncat,
		       &agreement, &split, rp.csplit);
	if (agreement <= count) continue;  /*no better than default */

	/*  sort it onto the list of surrogates */
	improve = (double)agreement / (me->primary)->count;  /* % agreement */
	ss = insert_split( &(me->surrogate), ncat, improve, rp.maxsur);
	if (ss !=0) {
	    ss->improve  = improve;
	    ss->var_num   = i;
	    ss->count     = 0;       /*corrected by nodesplit() */
	    if (rp.numcat[i]==0) {
		ss->spoint    = split;
		ss->csplit[0] = rp.csplit[0];
		}
	    else for (k=0; k<rp.numcat[i]; k++) ss->csplit[k] = rp.csplit[k];
	    }
	}
    }
