#SCCS %W% %G%
rpart.control <-
  function(minsplit=20, minbucket= round(minsplit/3), cp=.01,
	   maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10, 
	   surrogatestyle =0, maxdepth=30, ... ) {

	if (maxcompete<0) {
	    warning("The value of maxcompete supplied is <0; the value 0 was used instead")
	    maxcompete <-0
	    }
	if (any(xval<0)) {
	    warning("The value of xval supplied is <0; the value 0 was used instead")
	    xval <-0
	    }
	if (maxdepth > 30) stop("Maximum depth is 30")
	if (maxdepth < 1)  stop("Maximum depth must be at least 1")

	if (missing(minsplit) && !missing(minbucket)) minsplit <- minbucket*3

	# Because xval can be of length either 1 or n, and the C code
	#   refers to parameters by number, i.e., "opt[5]" in rpart.c, 
	#   the xval parameter should always be last on the list.
	list(minsplit=minsplit, minbucket=minbucket, cp=cp,
	     maxcompete=maxcompete, maxsurrogate=maxsurrogate,
	     usesurrogate=usesurrogate, 
	     surrogatestyle=surrogatestyle, maxdepth=maxdepth, xval=xval )
	}
