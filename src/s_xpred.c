/*
**  SCCS  @(#)s_xpred.c	1.8 02/08/98
** An S interface to "cross validated predictions"
**    99% of this routine is a copy of s_to_rp and rpart.c
*/
#include <setjmp.h>
#include <stdio.h>
#include "rpart.h"
#include "node.h"
#include "func_table.h"
#include "rpartS.h"
#include "rpartproto.h"

void s_xpred(long *sn, 	   long *nvarx,   long *ncat,    long *method, 
	     double *opt,  double *parms, long *xvals,   long *x_grp,
	     double *ymat, double *xmat,   long *missmat, double *predict,
	     long *ncp,    double *cp,    char **error)
    {
    int i,j,k;
    int maxcat;
    double temp;
    int n, nvar;
    int maxpri;
    struct node *xtree;
    int old_n;

    /*
    **  The opt string is in the order of control.rpart()
    **    minsplit, minbucket, cp, maxcomptete, maxsurrogate, usesurrogate,
    **    and xval
    */
    maxpri = opt[3] +1;

    /*
    ** Memory allocation errors from subroutines come back here
    */
    if (j=setjmp(errjump)) {
	*error = "Out of memory, cannot allocate needed structure";
	*sn = -1; 
	return;
	}

    /*
    ** initialize the splitting functions from the function table
    */
    if (*method <= NUM_METHODS) {
	i = *method -1;
	rp_init   = func_table[i].init_split;
	rp_choose = func_table[i].choose_split;
	rp_eval   = func_table[i].eval;
	rp_error  = func_table[i].error;
	rp.num_y  = func_table[i].num_y;
	}
    else {
	*error = "Invalid value for 'method'";
	*sn = -1; 
	return;
	}

    /*
    ** set some other parameters
    */
    n = *sn;
    nvar = *nvarx;
    rp.min_node =  opt[1];
    rp.min_split = opt[0];
    rp.complex   = opt[2];
    rp.nvar =  nvar;
    rp.numcat =  ncat;
    rp.maxpri = maxpri;
    if (maxpri <1) rp.maxpri =1;
    rp.maxsur = opt[4];
    rp.usesurrogate = opt[5];
    rp.n = n;
    rp.num_unique_cp = *ncp;

    /*
    ** create the "ragged array" pointers to the matrix
    **   x and missmat are in column major order
    **   y is in row major order
    */
    rp.xdata = (double **) ALLOC(nvar, sizeof(double *));
    if (rp.xdata==0) longjmp(errjump, 1);
    for (i=0; i<nvar; i++) {
	rp.xdata[i] = &(xmat[i*n]);
	}
    rp.ydata = (double **) ALLOC(n, sizeof(double *));
    if (rp.ydata==0) longjmp(errjump, 1);
    for (i=0; i<n; i++)  rp.ydata[i] = &(ymat[i*rp.num_y]);

    /*
    ** allocate some scratch
    */
    rp.tempvec = (int *)ALLOC(2*n, sizeof(int));
    rp.which   = rp.tempvec +n;
    rp.xtemp = (double *)ALLOC(n, sizeof(double));
    rp.ytemp = (double **)ALLOC(n, sizeof(double *));
    if (rp.tempvec==0 || rp.xtemp==0 || rp.ytemp==0) longjmp(errjump, 1);

    /*
    ** create a matrix of sort indices, one for each continuous variable
    **   This sort is "once and for all".  The result is stored on top
    **   of the 'missmat' array.
    ** I don't have to sort the categoricals.
    */
    rp.sorts  = (long**) ALLOC(nvar, sizeof(long *));
    maxcat=0;
    for (i=0; i<nvar; i++) {
	rp.sorts[i] = &(missmat[i*n]); 
	for (k=0; k<n; k++) {
	    if (rp.sorts[i][k]==1) {
		rp.tempvec[k] = -(k+1);
		rp.xdata[i][k] =0;     /* avoid weird numerics in S's NA */
	        }
	    else                   rp.tempvec[k] =  k;
	    }
	if (ncat[i]==0)  mysort(0, n-1, rp.xdata[i], rp.tempvec); 
	else if (ncat[i] > maxcat)  maxcat = ncat[i];
        for (k=0; k<n; k++) rp.sorts[i][k] = rp.tempvec[k];
	}

    /*
    ** And now the last of my scratch space
    */
    if (maxcat >0) {
	rp.csplit = (int *) ALLOC(3*maxcat, sizeof(int));
	if (rp.csplit==0) longjmp(errjump, 1);
	rp.left = rp.csplit + maxcat;
	rp.right= rp.left   + maxcat;
	}
    else rp.csplit = (int *)ALLOC(1, sizeof(int));

    (*rp_init)(n, rp.ydata, maxcat, error, parms, &rp.num_resp, 1);
    nodesize = sizeof(struct node) + (rp.num_resp-2)*sizeof(double);

    /*
    ** do the validations
    */
    old_n =n;
    for (i=0; i< *xvals; i++) {
	/*
	** mark the "leave out" data as ficticious node 0
	*/
	k=0;
	for (j=0; j<rp.n; j++) {
	    if (x_grp[j]==(i+1)) {
		rp.which[j] =0;
		}
	    else {
		rp.which[j] =1;
		rp.ytemp[k] = rp.ydata[j];
		k++;
		}
	    }

	/* rescale the cp */
	for (j=0; j<rp.num_unique_cp; j++) cp[j] *= (double)(k)/old_n;
	rp.alpha *= (double)k /old_n;
	old_n = k;

	/*
	** partition the new tree
	*/
	xtree = (struct node *) calloc(1, nodesize);
	xtree->num_obs = k;
	(*rp_init)(k,rp.ytemp, maxcat, error, parms, &temp, 2);
	(*rp_eval)(k, rp.ytemp, xtree->response_est, &(xtree->risk));
	xtree->complexity = xtree->risk;
	partition(1, xtree, &temp);
	fix_cp(xtree, xtree->complexity);
 
	/*
	** run the extra data down the new tree
	*/
	for (j=0; j<rp.n; j++) {
	    if (rp.which[j]==0) {
		rundown2(xtree, j, cp, (predict+ j* *ncp));
		}
	    }
	free_tree(xtree, 1);
	}
    }
