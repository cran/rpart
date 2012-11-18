/*
** prototypes for all of the rpart functions
**   This helps the ansi compiler do tight checking.
**
*/
struct node *branch(struct node *tree, int obs);

void bsplit(struct node *me, int n1, int n2);

void choose_surg(int n1, int n2, int *y,            double *x,     int *order, 
		 int ncat,       double *agreement, double *split, int *csplit,
		 double ltot,    double rtot,       double *adj);

void fix_cp(struct node *me, double parent_cp);

void free_tree(struct node *node,  int freenode);

void graycode_init0( int maxcat);
void graycode_init1( int numcat, int *count);
void graycode_init2( int numcat, int *count, double *val);
int  graycode(void);

struct split *insert_split(struct split **listhead, int ncat, 
			   double improve,          int max);

void make_cp_list(struct node *me, double parent, 
		  struct cptable *cptable_head);

struct cptable *make_cp_table(struct node *me, double parent, int nsplit);

void mysort(int start, int stop, double *x, int *cvec);

void nodesplit(struct node *me, int nodenum, int n1, int n2, 
	       int *nleft, int *nright);

int partition(int nodenum, struct node *splitnode, double *sumrisk,
	      int n1,      int n2);

int print_tree(struct node *me, int maxdepth);

void pyears2r(int   *sn,       int   *sny,    int   *sdoevent, 
	      double *sy,      double *wt,    int   *sodim,    int   *ofac, 
	      int   *odims,    double *socut, double *sodata,
	      double *pyears,  double *pn,    double *pcount, 
	      double *offtable);

SEXP rpart(SEXP ncat2,   SEXP method2,  SEXP opt2,
           SEXP parms2,  SEXP ymat2,    SEXP xmat2,
	   SEXP xvals2,  SEXP xgrp2,  	SEXP wt2,     
	   SEXP ny2,     SEXP cost2);

void rpart_callback0(int *nr);
void rpart_callback1(int n, double *y[], double *wt, double *z);
void rpart_callback2 (int n, int ncat, double *y[], double *wt, 
		     double *x, double *good);
void rpcountup(struct node *me, int *nnode, int *nsplit, int *ncat);

void rpmatrix(struct node *me, int *numcat,      double **dsplit,
	      int **isplit,    int **csplit,     double **dnode, 
	      int **inode,     int id);

void rundown(struct node *tree,  int obs,       double *cp, 
	     double *xpred,      double *xtemp);

void rundown2(struct node *tree, int obs, double *cp, double *xpred,
	      int nresp);

void surrogate(struct node *me, int n1, int n2);

SEXP xpred(SEXP ncat2,   SEXP method2,  SEXP opt2,
           SEXP parms2,  SEXP xvals2,   SEXP xgrp2,
	   SEXP ymat2,   SEXP xmat2,    SEXP wt2,
	   SEXP ny2,     SEXP cost2,    SEXP all2,
           SEXP cp2,     SEXP toprisk2, SEXP nresp2);

void xval(int n_xval,  struct cptable *cptable_head,  int *x_grp, 
	  int maxcat,  char **error,                  double * parms,
	  int *savesort);
