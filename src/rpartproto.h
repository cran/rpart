/*
**  SCCS @(#)rpartproto.h	1.1 02/08/98
** prototypes for all of the rpart functions
**   This helps the ansi compiler do tight checking.
**
*/

struct node *branch(struct node *tree, int obs);

void bsplit(struct node *me, int nodenum);

void choose_surg(int nodenum,    int *y,         double *x,     long *order, 
		 int ncat,       int *agreement, double *split, int *csplit);

void fix_cp(struct node *me, double parent_cp);

void free_tree(struct node *node,  int freenode);

void graycode_init0( int maxcat);
void graycode_init1( int numcat, int *count);
void graycode_init2( int numcat, int *count, double *val);
int  graycode();

struct split *insert_split(struct split **listhead, int ncat, 
			   double improve,          int max);

void make_cp_list(struct node *me, double parent, 
		  struct cptable *cptable_head);

struct cptable *make_cp_table(struct node *me, double parent, int nsplit);

void mysort(int start, int stop, double *x, int *cvec);

void nodesplit(struct node *me, int nodenum);

int partition(int nodenum, struct node *splitnode, double *sumrisk);

void pred_rpart(long *dimx,	long *nnode, 	long *nsplit, 	long *dimc, 
		long *nnum,  	long *nodes2,   long *vnum,     double *split2,
		long *csplit2,  long *usesur,   double *xdata2, 
		long *xmiss2,   long *where);

int rpart(int n,         int nvarx,      long *ncat,     int method, 
          int mnode,     int msplit,     int  maxpri,    int maxsur,
	  int usesur,    double *parms,  double *ymat,   double *xmat,  
          long *missmat, double complex, struct cptable *cptable,
	  struct node **tree,            char **error,   int *which,
	  int xvals,     long *x_grp);

void rpcountup(struct node *me, long *nnode, long *nsplit, int *ncat);

void rplabel(long *nsplit,   long *index,   double *splits, 
             long *ncat,     long *csplit,  char   **cutleft, char **cutright);

void rpmatrix(struct node *me,  long *nodecount,   long *splitcount, 
	      long *catcount,   long *numcat,      double **dsplit,
	      long **isplit,    long **csplit,     double **dnode, 
	      long **inode,     int id);

void rundown(struct node *tree,  int obs,     double *cp, 
	     double *xpred,      double *xtemp);

void rundown2(struct node *tree, int obs, double *cp, double *xpred);

void s_to_rp(long *n, 	  long *nvarx, 	 long *ncat, 	long *method, 
	     double *opt, double *parms, long *xvals,   long *x_grp,
	     double *y,   double *xmat,   long *missmat, char **error);

void s_to_rp2(long *n,         long *nsplit,    long *nnode,     long *ncat, 
	      long *numcat,    long *maxcat,    long *xvals,     long *which, 
	      double *cptable, double *dsplit,  long *isplit,    long *csplit,
	      double *dnode,   long *inode);

void s_xpred(long *sn, 	   long *nvarx,   long *ncat,    long *method, 
	     double *opt,  double *parms, long *xvals,   long *x_grp,
	     double *ymat, double *xmat,   long *missmat, double *predict,
	     long *ncp,    double *cp,    char **error);

void surrogate(struct node *me, int nodenum);

void xval(int n_xval,  struct cptable *cptable_head,  long *x_grp, 
	  int maxcat,  char **error,                  double * parms);
