/* SCCS @(#)rpart.h	1.5 12/13/99 */
/*
** commom variables for the rpart routine
*/
#ifndef FLOAT
#define FLOAT double    /*so I can easily change 'x' to double at some later
			 date, with all the consequences thereof.  Also see
		         node.h  */
#endif
#define LEFT  (-1)     /*used for the variable "extra" in nodes */
#define RIGHT  1
#define MISSING 0
#define CPLIST_EPS .05  /* Any complexity parameters within 10% of each other
			**   are considered tied, i.e., (x-y)/(x+y) < .05 .
			*/

/* As a sop to S, I need to keep the total number of external symbols
**  somewhat smaller.  So, pack most of them all into a structure.
*/
struct {
    double complex;
    double alpha;
    double **ydata;
    FLOAT  **xdata;
    FLOAT  *xtemp;
    double *wt;
    double **ytemp;
    double *wtemp;          /* temp vector of weights */
    double *lwt;
    double *rwt;            /*scratch double vectors, of length ncat */
    int    *numcat;        /* variable type: 0=cont, 1+  =#categories */
    int   **sorts;             /* allocated on the fly */
    int    n;              /* total number of subjects  */
    int    num_y;          /* number of y variables */
    int    nvar;           /* number of predictors */
    int    maxpri;
    int    maxsur;   /* max # of primary or surrogate splits to use */
    int    usesurrogate;
    int    num_unique_cp;
    int    min_node;       /* minimum size for any terminal node */
    int    min_split;      /*minimum size before we attempt a split */
    int    num_resp;       /*length of the response vector */
    int    sur_agree;      /*0 =  my style, 1=CART style */
    int    *tempvec;       /*to be allocated by the mainline, of length n */
    int    *which;
    int    *csplit;
    int    *left;
    int    *right;
    }  rp;

struct cptable *cptable_tail;
int  (*rp_init)();    /*called to initialize a splitting function */
void (*rp_choose)();  /*set to the splitting function */
void (*rp_eval)() ;   /*set to the evaluation routine */
double (*rp_error)();     /*set to the prediction error routine */

/*
** The user inputs his complexity parameter as a percentage. and the
**   printout is also scaled in this way.  The book and the computations all
**   have an easier time with absolute cp.  So complex_p = what the user
**   typed and alpha = complex_p * (risk of top node) = what is used
**   internally.
**
**
** Categorical variables must be coded as 1,2,3, ..., and there may be
**  missing categories.  The upper limit is determined on the fly.
*/
