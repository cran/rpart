/* SCCS @(#)rpmatrix.c	1.5 02/08/98 */
/*
**  For S's usage, convert the linked list data into matrix form
*/
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

void rpmatrix(struct node *me,  long *nodecount,   long *splitcount, 
	      long *catcount,   long *numcat,      double **dsplit,
	      long **isplit,    long **csplit,     double **dnode, 
	      long **inode,     int id)
    {
    /*
    ** dsplit  0: improvement
    **         1: split point if continuous; index into csplit if not
    ** isplit  0: variable #
    **         1: count
    **         2: if continuous: direction -1=left, 1=right
    **            if categorical: # of categories
    ** csplit[i]: -1=left, 0=missing category, 1=right
    ** dnode   0: risk
    **         1: complexity threshold
    **         2-?: response estimate
    ** inode   0: node number
    **         1: index of the first primary, in the split list
    **         2: #primary    ==0 if this is a terminal node
    **         3: #surrogates
    **         4: # observations
    **         5: # obs for which this is the final resting place
    */

    int i,j, k;
    struct split *spl;
    static int scnt, ncnt, ccnt;
    static double cp_scale;

    if (id==1) {
	cp_scale = 1/ me->risk;
	}
    scnt = *splitcount;
    ncnt = *nodecount;
    ccnt = *catcount;
    dnode[0][ncnt] = me->risk;
    dnode[1][ncnt] = me->complexity * cp_scale;
    for (i=0; i<rp.num_resp; i++) dnode[2+i][ncnt] = me->response_est[i];
    inode[0][ncnt] = id;
    inode[4][ncnt] = me->num_obs;

    if (me->complexity <=rp.alpha || me->leftson==0) { /*no kids */
	inode[1][ncnt] = 0;
	inode[2][ncnt] = 0;
	inode[3][ncnt] = 0;
	inode[5][ncnt] = me->num_obs;
	*nodecount = ncnt+1;
	}
    else {
	inode[1][ncnt] = scnt +1;    /*S has 1 based, not 0 based subscripts */

	i=0;
	for (spl = me->primary; spl!=0; spl = spl->nextsplit) {
	    i++;
	    j = spl->var_num;
	    dsplit[0][scnt] = spl->improve;
	    if (numcat[j] ==0) {
		dsplit[1][scnt] = spl->spoint;
		isplit[2][scnt] = spl->csplit[0];
		}
	    else {
		dsplit[1][scnt] = ccnt+1;
		isplit[2][scnt] = numcat[j];
		for (k=0; k<numcat[j]; k++) csplit[k][ccnt] = spl->csplit[k];
		ccnt++;
		}
	    isplit[0][scnt] = j +1;      /* use "1" based subscripts */
	    isplit[1][scnt] = spl->count;
	    scnt++;
	    }
	inode[2][ncnt] = i;

	i=0;
	for (spl=me->surrogate; spl!=0; spl = spl->nextsplit) {
	    i++;
	    j = spl->var_num;
	    dsplit[0][scnt] = spl->improve;
	    if (numcat[j] ==0) {
		dsplit[1][scnt] = spl->spoint;
		isplit[2][scnt] = spl->csplit[0];
		}
	    else {
		dsplit[1][scnt] = ccnt+1;
		isplit[2][scnt] = numcat[j];
		for (k=0; k<numcat[j]; k++) csplit[k][ccnt] = spl->csplit[k];
		ccnt++;
		}
	    isplit[0][scnt] = j +1;
	    isplit[1][scnt] = spl->count;
	    scnt++;
	    }
	inode[3][ncnt] = i;
	inode[5][ncnt] = me->num_obs -
			   ((me->leftson)->num_obs + (me->rightson)->num_obs);

	ncnt++;
	*nodecount = ncnt;
	*splitcount= scnt;
	*catcount  = ccnt;

	rpmatrix(me->leftson, nodecount, splitcount, catcount, numcat,
		    dsplit, isplit, csplit, dnode, inode, 2*id);
	rpmatrix(me->rightson,nodecount, splitcount, catcount, numcat,
		    dsplit, isplit, csplit, dnode, inode, 2*id +1);
	}
    }
