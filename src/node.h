/* SCCS @(#)node.h	1.3 12/13/99 */
/*
** definition of a node in the tree
*/
#ifndef FLOAT
#define FLOAT float   /*see comments in rpart.h */
#endif

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif

struct split {
    double improve;
    double adj;           /* for surrogates only, adjusted agreement */
    FLOAT spoint;    /*only used if it is continuous */
    struct split *nextsplit;
    int var_num;
    int count;
    int csplit[2];     /*the actual length will be longer for a categorical */
    };                 /*   predictor with >2 levels */

struct node {
    double  risk;       /*risk for the node */
    double complexity;  /* complexity at which it will collapse */
    double sum_wt;      /* sum of the weights for the node  */
    struct split *primary;
    struct split *surrogate;
    struct node *rightson;
    struct node *leftson;
    int num_obs;
    int lastsurrogate;
    double response_est[2];  /*actual length depends on splitting rule */
    };
EXTERN int nodesize;

struct cptable {
    double cp;
    double risk;
    double xrisk;
    double xstd;
    int nsplit;
    struct cptable *forward;
    struct cptable *back;
    }  ;

/**************************************************************************
**
**  Split:
**      variable number of the split; 0 = no more surrogates (or primaries)
**
**      split point: the actual split point for a continuous
**
**      improve:  For primary splits, the iprovement index returned by the
**                 bsplit routine.  This is the measure that determines the
**                 winning split.
**                For surrogate splits, this holds the error rate, i.e., the
**                 % incorrect guesses of the primary by using this surrogate.
**
**      count: The number of observations split using this variable.  For the
**             first primary, this will = the number of non-missing values.
**             For surrogates, it will be the number missing in the primary
**             and all earlier surrogates but not missing on this one.  (For
**             all primaries but the first, the number is theoretical).
**
**	adj:  Let "maj" be the %agreement for going with the majority,
**                and "agree" the %agreement for this surrogate.  The
**                adjusted value is (agree - maj)/(1-maj); the amount of
**                the potential improvement actually realized.  The denominator
**                for both percents depends on the sur_agree option.
**
**      csplit[0]:   For a continuous variable, we also need to know the
**                    direction of the split.  We use this "extra" variable
**                    as 1: <x to the left, -1: <x to the right.
**
**      csplit[]:    For a categorical, the labels are LEFT, RIGHT, and
**                    0=missing.  (Even if a particular category is not empty,
**                    there may be no subjects from that category present
**                    at a particular split further down the tree).
**
**
**  Node:
**      num_obs: Number of observations in the node.
**
**      response_est: From the eval routine.  Estimate of the response, if
**                      this node were terminal.
**
**      risk: From the eval routine. Estimate of risk, if this node were
**                      terminal.
**
**      complexity: On the way down, it holds equation 5.18.  On the way up
**              it holds a provisional C.P. (The actual C.P. for each node
**              will be the minimum of this number and the provisional C.P.
**              of all nodes above it.  One more pass downward can establish
**              the proper C.P.).
**
**      lastsurrogate: Which direction to send obs for which the primary and
**              all the surrogates are missing.  (The child with the greatest
**		sum of weights).
*/
