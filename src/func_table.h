/* SCCS @(#)func_table.h	1.3  01/25/97 */
/*
** The table of implimented splitting functions
**
**  init_split   - Will be called before a tree is started.  May do very
**                  little, but empirical Bayes like methods will need to set
**                  some global variables.
**  choose_split - function to find the best split
**  eval         - function to calculate the response estimate and risk
**  error        - Function that returns the prediction error.
**  num_y        - Number of columns needed to represent y (usually 1)
*/

extern void anova(), anovass();
extern int anovainit();
extern double anovapred();
extern int poissoninit();
extern void poisson(), poissondev();
extern double poissonpred();
extern int giniinit();
extern void gini(), ginidev();
extern double ginipred();

static struct {
	int  (*init_split)();
	void (*choose_split)();
	void (*eval)();
	double (*error)();
	int num_y;
	  } func_table [] =
		 {{anovainit,   anova,   anovass,    anovapred, 1},
		  {poissoninit, poisson, poissondev, poissonpred, 2},
		  {giniinit, gini, ginidev, ginipred, 1},
		 };

#define NUM_METHODS 3  /*size of the above structure */
