/* SCCS @(#)anova.c	1.4  02/08/98 */
/*
** The four routines for anova splitting
*/
#include "rpartS.h"
#define LEFT  -1     /*used for the variable "extra" in nodes */
#define RIGHT  1
#define MISSING 0

static double *mean;
static    int *countn, *order, *order2;

int anovainit(int n,        double *y[],  int maxcat, char **error, 
	      double *parm, int *size,    int who)
    {
    if (who==1 && maxcat >0) {
	countn = (int *)ALLOC(3*maxcat, sizeof(int));
	order  = countn + maxcat;
	order2 = order + maxcat;
	mean  = (double *)ALLOC(maxcat, sizeof(double));
	if (countn==0 || mean==0) {
	    *error = "Could not allocate memory in anovainit";
	    return(1);
	    }
	}
    *size =1;
    return(0);
    }
/*
** The anova evaluation function.  Return the mean and the ss.
*/
void anovass(int n, double *y[], double *value, double *risk)
    {
     int i;
     double temp;
    double mean, ss;

    temp =0;
    for (i=0; i<n; i++)
	temp += *y[i];

    mean = temp/n;

    ss =0;
    for (i=0; i<n; i++) {
	temp = *y[i] - mean;
	ss += temp*temp;
	}

    *value = mean;
    *risk = ss;
    }

/*
** The anova splitting function.  Find that split point in x such that
**  the sum of squares of y within the two groups is decreased as much
**  as possible.  It is not necessary to actually calculate the SS, the
**  improvement involves only means in the two groups.
*/
void anova(int n,    double *y[],     double *x,     int nclass, 
	   int edge, double *improve, double *split, int *csplit, double myrisk)
    {
    int i,j;
    double temp;
    double left_sum, right_sum;
    int left_n, right_n;
    double grandmean, best;
    int direction;
    int where;
    int ncat;

    /*
    ** Compute the grand mean, which will be subtracted from all of the
    **  data elements "on the fly".  This makes the hand calculator formula
    **  numerically stable.
    ** Also get the total n
    */
    right_n =n;
    grandmean =0;
    for (i=0; i<n; i++) grandmean += *y[i];
    grandmean /= right_n;

    if (nclass==0) {
	right_sum =0; left_sum=0;
	left_n =0;
	best    = 0;
	for (i=0; right_n>edge; i++) {
	    left_n++;  right_n--;
	    temp = *y[i] - grandmean;
	    left_sum  +=temp;
	    right_sum -=temp;
	    if (x[i+1] !=x[i] &&  left_n>=edge) {
		temp = left_sum*left_sum/left_n  +
			    right_sum*right_sum/right_n;
		if (temp > best) {
		    best=temp;
		    where =i;
		    if (left_sum < right_sum) direction = LEFT;
				      else    direction = RIGHT;
		    }
		}
	    }

	*improve =  best/ myrisk;
	if (best>0) {   /* found something */
	    csplit[0] = direction;
	    *split = (x[where] + x[where+1]) /2;
	    }
	}

    else {
	for (i=0; i<nclass; i++) {
	    mean[i] =0;
	    countn[i]=0;
	    }

	for (i=0; i<n; i++) {
	    j = x[i] -1;
	    countn[j]++;
	    mean[j] += *y[i] - grandmean;
	    }

	/*
	** Rank the means  - each is scored as the number of others that it
	**  is smaller than.  Ignore the categories which had no representatives.
	*/
	for (i=0; i<nclass; i++) {
	    order[i]=0;
	    if (countn[i]>0) {
		mean[i] /= countn[i];
		for (j=i-1; j>=0; j--) {
		    if (countn[j] >0) {
			if (mean[i]>mean[j])  order[j]++;
			    else              order[i]++;
			}
		    }
		}
	    }
	/*
	** order2 will point to the largest, second largest, etc
	*/
	ncat =0;
	for (i=0; i<nclass; i++) {
	    if (countn[i]>0) {
		ncat++;
		order2[order[i]]=i;
		}
	    }

	/*
	** Now find the split that we want
	*/
	left_n =0;  right_n = n;
	left_sum=0; right_sum =0;
	best =0;
	where =0;
	for (i=0; i<ncat-1; i++){
	    j = order2[i];
	    left_n += countn[j];
	    right_n-= countn[j];
	    left_sum += countn[j]*mean[j];
	    right_sum-= countn[j]*mean[j];
	    if (left_n>=edge  &&  right_n>=edge) {
		temp = left_sum*left_sum/left_n  + right_sum*right_sum/right_n;
		if (temp > best) {
		    best=temp;
		    where =i;
		    if (left_sum > right_sum) direction = RIGHT;
				      else    direction = LEFT;
		    }
		}
	    }

	*improve = best/ myrisk;      /* % improvement */

	/* (if best=0, csplit will never be looked at by the calling routine) */
	for (i=0; i<nclass; i++)      csplit[i]=0;
	for (i=0; i<=where; i++)  csplit[order2[i]] = direction ;
	for (   ; i<ncat; i++)  csplit[order2[i]] =  -direction;
	}
    }
