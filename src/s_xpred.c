/*
** An S interface to "cross validated predictions"
**    99% of this routine is a copy of s_to_rp and rpart.c and xval.c
** Note that nresp is the length of the response vector that is desired
**   and will be returned, not necessarily the actual length.  It is common
**   to only request the first element of a long response vector.
*/
#include "rpart.h"
#include "node.h"
#include "func_table.h"
#include "rpartproto.h"

void s_xpred(int *sn,	   int *nvarx,    int *ncat,    int *method,
	     double *opt,  double *parms, int *xvals,   int *x_grp,
	     double *ymat, double *xmat,   double *predict,
	     int *ncp,     double *cp,    char **error,  double *wt,
	     int *ny,      double *cost,  int *nresp2)
{
    int i,j,k;
    int ii;
    int last;
    int maxcat;
    int nresp;
    double temp;
    int n, nvar;
    int maxpri;
    struct node *xtree;
    double old_wt, total_wt;
    int *sj;
    int *si, *savesort;
    int xgroup;

    /*
    ** initialize the splitting functions from the function table
    */
    if (*method <= NUM_METHODS) {
	i = *method -1;
	rp_init   = func_table[i].init_split;
	rp_choose = func_table[i].choose_split;
	rp_eval   = func_table[i].eval;
	rp_error  = func_table[i].error;
	rp.num_y  = *ny;
    }
    else {
	*error = _("Invalid value for 'method'");
	*sn = -1;
	return;
    }

    /*
    ** Set other parameters
    **
    **  The opt string is in the order of control.rpart()
    **    minsplit, minbucket, cp, maxcomptete, maxsurrogate, usesurrogate,
    **    and xval
    */
    n = *sn;
    nresp = *nresp2;
    nvar = *nvarx;
    rp.nvar =  nvar;
    rp.numcat = ncat;
    rp.n = n;
    rp.num_unique_cp = *ncp;
    rp.wt = wt;

    rp.min_node =  (int)opt[1];
    rp.min_split = (int)opt[0];
    rp.complexity= opt[2];
    maxpri       = (int)opt[3] +1;
    rp.maxpri = maxpri;
    if (maxpri <1) rp.maxpri =1;
    rp.maxsur = (int)opt[4];
    rp.usesurrogate = (int)opt[5];
    rp.sur_agree = (int)opt[6];
    rp.maxnode  = (int)pow((double)2.0, opt[7]) -1;
    rp.vcost    = cost;

    /*
    ** create the "ragged array" pointers to the matrix
    **   x and missmat are in column major order
    **   y is in row major order
    */
    rp.xdata = (double **) ALLOC(nvar, sizeof(double *));
    for (i=0; i<nvar; i++) {
	rp.xdata[i] = &(xmat[i*n]);
    }
    rp.ydata = (double **) ALLOC(n, sizeof(double *));
    for (i=0; i<n; i++)  rp.ydata[i] = &(ymat[i*rp.num_y]);

    /*
    ** allocate some scratch
    */
    rp.tempvec = (int *)ALLOC(2*n, sizeof(int));
    rp.which   = rp.tempvec +n;
    rp.xtemp = (double *)ALLOC(n, sizeof(double));
    rp.ytemp = (double **)ALLOC(n, sizeof(double *));
    rp.wtemp = (double *)ALLOC(n, sizeof(double));

    /*
    ** create a matrix of sort indices, one for each continuous variable
    **   This sort is "once and for all".  The result is stored on top
    **   of the 'missmat' array.
    ** I don't have to sort the categoricals.
    */
    rp.sorts  = (int**) ALLOC(nvar, sizeof(int *));
    rp.sorts[0]=(int *) ALLOC(nvar*n, sizeof(int));
    maxcat=0;
    for (i=0; i<nvar; i++) {
	rp.sorts[i] = rp.sorts[0] + i*n;
	for (k=0; k<n; k++) {
	    if (!R_FINITE(rp.xdata[i][k])) {
		rp.tempvec[k] = -(k+1);
		rp.xtemp[k] =0;     /* avoid weird numerics in S's NA */
	    }
	    else  {
		rp.tempvec[k] =  k;
		rp.xtemp[k] = rp.xdata[i][k];
	    }
	}
	if (ncat[i]==0) {
	    /* We want to leave the x matrix alone, and only get indices */
	    mysort(0, n-1, rp.xtemp, rp.tempvec);
	}
	else if (ncat[i] > maxcat)  maxcat = ncat[i];
	for (k=0; k<n; k++) rp.sorts[i][k] = rp.tempvec[k];
    }

    /*
    ** save away a copy of the rp.sorts
    */
    savesort = (int*) ALLOC(n*nvar, sizeof(int));
    si = savesort;
    sj = rp.sorts[0];
    for (i=0; i < n *  rp.nvar; i++) *si++ = *sj++;

    /*
    ** And now the last of my scratch space
    */
    if (maxcat >0) {
	rp.csplit = (int *) ALLOC(3*maxcat, sizeof(int));
	rp.left = rp.csplit + maxcat;
	rp.right= rp.left   + maxcat;
	rp.lwt    = (double *) ALLOC(2*maxcat, sizeof(double));
	rp.rwt  = rp.lwt    + maxcat;
    }
    else rp.csplit = (int *)ALLOC(1, sizeof(int));

    (*rp_init)(n, rp.ydata, maxcat, error, parms, &rp.num_resp, 1, rp.wt);
    nodesize = sizeof(struct node) + (rp.num_resp-2)*sizeof(double);

    /*
    ** I need the risk of the full tree, to scale alpha
    */
    xtree = (struct node *) ALLOC(1, nodesize);
    memset(xtree, 0, sizeof(struct node));
    (*rp_eval)(n, rp.ydata, xtree->response_est, &(xtree->risk), rp.wt);
    rp.alpha = rp.complexity * (xtree)->risk;

    free_tree(xtree, 0);

    /*
    ** do the validations
    */
    total_wt = 0;
    for (i = 0; i < rp.n; i++) total_wt += rp.wt[i];
    old_wt = total_wt;

    k = 0; /* -Wall */
    for (xgroup=0; xgroup<*xvals; xgroup++) {
	/*
	** restore rp.sorts, with the data for this run at the top
	**   this requires one pass per variable
	*/
	for (j=0; j<rp.nvar; j++) {
	    k=0;
	    for (i=0; i<rp.n; i++) {
		ii = savesort[j*n +i];  /* walk through the variable in order*/
		if (ii<0) ii = -(1+ii);  /* missings move too */
		if (x_grp[ii]!=(xgroup+1)) {
		    /*
		    ** this obs is left in --
		    **  copy to the front half of rp.sorts
		    */
		    rp.sorts[j][k] = savesort[j*n + i];
		    k++;
		}
	    }
	}

	/*
	**  Fix up the y vector, and save a list of "left out" obs
	**   in the tail, unused end of rp.sorts[0][i];
	*/
	last=k;
	k=0;
	temp =0;
	for (i=0; i<n; i++) {
	    if (x_grp[i] == (xgroup +1)) {
		rp.sorts[0][last] = i;
		last++;
	    }
	    else {
		rp.ytemp[k] = rp.ydata[i];
		rp.wtemp[k] = rp.wt[i];
		temp += rp.wt[i];
		k++;
	    }
	}

	/* at this point k = #obs in the xval group */
	/* rescale the cp */
	for (j=0; j<rp.num_unique_cp; j++) cp[j] *= temp/old_wt;
	rp.alpha *= temp/old_wt;
	old_wt = temp;

	/*
	** partition the new tree
	*/
	xtree = (struct node *) CALLOC(1, nodesize);
	xtree->num_obs = k;
	(*rp_init)(k,rp.ytemp, maxcat, error, parms, &temp, 2, rp.wtemp);
	(*rp_eval)(k, rp.ytemp, xtree->response_est, &(xtree->risk),
		   rp.wtemp);
	xtree->complexity = xtree->risk;
	partition(1, xtree, &temp, 0, k);
	fix_cp(xtree, xtree->complexity);

	/* if (xgroup==0) print_tree(xtree, 1, 0, 0, 0); debug line */
	/*
	** run the extra data down the new tree
	*/
	for (i = k; i < rp.n; i++) {
	    j = rp.sorts[0][i];
	    rundown2(xtree, j, cp, (predict+ j* (*ncp) * nresp), nresp);
	}

	free_tree(xtree, 1);  // Calloc-ed inside loop
    }
}
