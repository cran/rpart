#SCCS @(#)rpart.control.s	1.4 01/31/97
rpart.control <-
  function(minsplit=20, minbucket= round(minsplit/3), cp=.01,
	   maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10,
	   surrogatestyle =0,  ... ) {

	if (maxcompete<0) {
	    warning("The value of maxcompete supplied is <0; the value 0 was used instead")
	    maxcompete <-0
	    }
	if (any(xval<0)) {
	    warning("The value of xval supplied is <0; the value 0 was used instead")
	    xval <-0
	    }

	if (missing(minsplit) && !missing(minbucket)) minsplit <- minbucket*3

	list(minsplit=minsplit, minbucket=minbucket, cp=cp,
	     maxcompete=maxcompete, maxsurrogate=maxsurrogate,
	     usesurrogate=usesurrogate,
	     surrogatestyle=surrogatestyle, xval=xval )
	}
