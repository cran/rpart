/* SCCS @(#)insert_split.c	1.3 02/08/98
/*
** sort a new split into a linked list, based on its "improvement"
**
**  allocates new memory as needed
**   returns 0 if the new element isn't good enough,
**   the address of the new element otherwise
*/
#include "node.h"
#include "rpart.h"
#include "rpartproto.h"

struct split *insert_split(struct split **listhead, int ncat, 
			   double improve,          int max)
    {
    int nlist;
    struct split *s1, *s2, *s3, *s4;

    if (ncat==0) ncat++;     /*number of ints needed in the structure */
    if (*listhead ==0) {
	/* first call to a new list */
	s3 = (struct split *)calloc(1, sizeof(struct split)+
						(ncat-1)*sizeof(int));
	if (s3==0) longjmp(errjump, 2);   /*no memory */
	s3->nextsplit =0;
	*listhead = s3;
	return(s3);
	}

    if (max <2) {
	/* user asked for only 1 to be retained! */
	s3 = *listhead;
	if (improve <= s3->improve) return(0);
	if (ncat >1) {
	    free(s3);
	    s3 = (struct split *)calloc(1, sizeof(struct split)+
						(ncat-1)*sizeof(int));
	    if (s3==0) longjmp(errjump, 2);   /*no memory */
	    s3->nextsplit =0;
	    *listhead = s3;
	    }
	return(s3);
	}

    /*set up --- nlist = length of list, s4=last element, s3=next to last */
    nlist=1;
    for (s4= *listhead; s4->nextsplit!=0; s4=s4->nextsplit) {
	s3 = s4;
	nlist++;
	}

    /* now set up so that the "to be added" is between s1 and s2 */
    s1 = *listhead;
    for (s2= *listhead; s2 !=0; s2= s2->nextsplit) {
	if (improve > s2->improve) break;
	s1 = s2;
	}

    if (nlist==max) {
	if (s2==0)  return(0);        /* not good enough */
	if (ncat >1) {
	    free(s4);              /*get new memory-- this chunk may be too small */
	    s4 = (struct split *)calloc(1, sizeof(struct split) +
						     (ncat-2)*sizeof(int));
	    if (s4==0) longjmp(errjump,2);
	    }
	if (s1==s3)   s4->nextsplit =0;
	else         {s3->nextsplit =0;  s4->nextsplit =s2;}
	}
    else {
	s4 = (struct split *)calloc(1, sizeof(struct split) +
							(ncat-2)*sizeof(int));
	if (s4==0) longjmp(errjump, 2);   /*no memory */
	s4->nextsplit = s2;
	}
    if (s2== *listhead) *listhead     = s4;
    else                s1->nextsplit = s4;
    return(s4);
    }
