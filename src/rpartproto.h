/*
**  SCCS  @(#)rpartproto.h	1.2 12/13/99
** prototypes for all of the rpart functions
**   This helps the ansi compiler do tight checking.
**
*/

struct node *branch(struct node *tree, int obs);

void bsplit(struct node *me, int nodenum);

void choose_surg(int nodenum,    int *y,            FLOAT *x,     int *order, 
		 int ncat,       double *agreement, FLOAT *split, int *csplit,
		 double ltot,    double rtot,       double *adj);

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

void mysort(int start, int stop, FLOAT *x, int *cvec);

void nodesplit(struct node *me, int nodenum);

int partition(int nodenum, struct node *splitnode, double *sumrisk);

void pred_rpart(int *dimx,	int *nnode, 	int *nsplit, 	int *dimc, 
		int *nnum,  	int *nodes2,   int *vnum,     double *split2,
		int *csplit2,  int *usesur,   double *xdata2, 
		int *xmiss2,   int *where);

int rpart(int n,         int nvarx,      int *ncat,     int method, 
          int mnode,     int msplit,     int  maxpri,    int maxsur,
	  int usesur,    double *parms,  double *ymat,   FLOAT *xmat,  
          int *missmat, double complex, struct cptable *cptable,
	  struct node **tree,            char **error,   int *which,
	  int xvals,     int *x_grp,    double *wt,     int surragree);

void rpcountup(struct node *me, int *nnode, int *nsplit, int *ncat);

void rplabel(int *nsplit,   int *index,   double *splits, 
             int *ncat,     int *csplit,  char   **cutleft, char **cutright);

void rpmatrix(struct node *me,  int *nodecount,   int *splitcount, 
	      int *catcount,   int *numcat,      double **dsplit,
	      int **isplit,    int **csplit,     double **dnode, 
	      int **inode,     int id);

void rundown(struct node *tree,  int obs,     double *cp, 
	     double *xpred,      double *xtemp);

void rundown2(struct node *tree, int obs, double *cp, double *xpred);

void s_to_rp(int *n, 	  int *nvarx, 	 int *ncat, 	int *method, 
	     double *opt, double *parms, int *xvals,   int *x_grp,
	     double *y,   FLOAT  *xmat,  int *missmat, char **error,
	     double *wt);

void s_to_rp2(int *n,         int *nsplit,    int *nnode,     int *ncat, 
	      int *numcat,    int *maxcat,    int *xvals,     int *which, 
	      double *cptable, double *dsplit,  int *isplit,    int *csplit,
	      double *dnode,   int *inode);

void s_xpred(int *sn, 	   int *nvarx,   int *ncat,    int *method, 
	     double *opt,  double *parms, int *xvals,   int *x_grp,
	     double *ymat, FLOAT  *xmat,  int *missmat, double *predict,
	     int *ncp,    double *cp,    char **error,  double *wt);

void surrogate(struct node *me, int nodenum);

void xval(int n_xval,  struct cptable *cptable_head,  int *x_grp, 
	  int maxcat,  char **error,                  double * parms);
