#SCCS  @(#)printcp.s	1.5 09/03/97
# print out the cptable, along with some summary of the tree
printcp <- function(x, digits=.Options$digits-2) {
    if (!inherits(x, 'rpart')) stop ("Must be an rpart x")
    cat(switch(x$method,anova = "\nRegression tree:\n" ,
			class = "\nClassification tree:\n" ,
			poisson="\nRates regression tree:\n",
			exp = "\nSurvival regression tree:\n")
			)

    if(!is.null(cl <- x$call)) {
	dput(cl)
	cat("\n")
      }
    frame <- x$frame
    leaves <- frame$var == "<leaf>"
    used <- unique(frame$var[!leaves])

    if(!is.null(used)) {
		cat("Variables actually used in tree construction:\n")
		print(sort(as.character(used)), quote=F)
		cat("\n")
	}


    cat("Root node error: ", format(frame$dev[1], digits=digits), '/', 
			frame$n[1], ' = ',
		         format(frame$dev[1]/frame$n[1], digits=digits), 
 			'\n\n', sep='')

    print (x$cptable, digits=digits)
    invisible(x$cptable)
    }

