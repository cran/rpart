/* SCCS @(#)make_cp_list.c	1.3 12/13/99 */
/*
** This routine creates the list of unique complexity parameters.
** The list is maintained in sorted order.  If two parameters are within
** "cplist_epsilon" of one another, then only the larger of them is
** retained.
**   I want the list sorted with the largest cp at the top of the list, since
** that is the order that the CP table will be printed in.
**
**    After the partition routine is done, each node is labeled with the
** complexity parameter appropriate if that node were the top of the tree.
** However, if there is a more vunerable node further up, the node in
** question will actually have the smaller complexity parameter; it will
** be removed when its parent collapses. So this routine also adjusts each
** C.P. to = minimum(my C.P., parent's C.P.).
**
**   This routine is called at the top level by rpart, after rpart has
** initialized the first member of the linked cp-list, set its number of
** splits to zero, and its risk to that for no splits at all.  This routine
** allocates and links in the rest of the cp-list.  The make_cp_table
** routine then fills in the rest of the variables in the list.
**
**  node *me;         pointer to my node structure  
**  double parent;    complexity of my parent node 
**
**   When it comes time to cross-validate, we fill in xrisk and xstd
*/
#include <math.h>
#include "rpart.h"
#include "node.h"
#include "rpartS.h"
#include "rpartproto.h"

void make_cp_list(struct node *me, double parent, struct cptable *cptable_head)
    {
    double me_cp;
    struct cptable *cplist, *cptemp;

    if (me->complexity > parent) me->complexity = parent;
    me_cp = me->complexity;
    if (me_cp < rp.alpha) me_cp = rp.alpha;    /*table should go no lower */
    if (me->leftson != 0) {
	 make_cp_list(me->leftson, me_cp, cptable_head);
	 make_cp_list(me->rightson,me_cp, cptable_head);
	 }

    if (me_cp < parent) {  /*if not, then it can't be unique */
	for (cplist= cptable_head; cplist !=0; cplist= cplist->forward) {
	    /* am I tied? */
	    if (me_cp==0 && cplist->cp==0) return;
	    if ((fabs(cplist->cp -me_cp)/ (cplist->cp + me_cp)) < CPLIST_EPS){
		if (me_cp > cplist->cp) cplist->cp = me_cp;
		return;
		}


	    /* should I be inserted? */
	    if (me_cp > cplist->cp) break;
	    cptemp = cplist;
	    }

	/* insert new stuff after cptemp */
	cplist = (struct cptable *) CALLOC(1, sizeof(struct cptable));
	cplist->cp = me_cp;
	cplist->xrisk = 0;
	cplist->xstd  =0;
	cplist->back = cptemp;
	cplist->forward = cptemp->forward;
	if (cptemp->forward!=0) (cptemp->forward)->back = cplist;
	  else  cptable_tail = cplist;
	cptemp->forward = cplist;
	rp.num_unique_cp++;
	return;
	}
    }