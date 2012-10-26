\name{NEWS}
\title{NEWS file for the rpart package}
\section{Changes in version 4.1-1, 2012 July 23}{
  \itemize{
    \item Merged in a set of Splus code changes that had accumulated
    at Mayo over the course of a decade. The primary one is
    a change in how indexing is done in the underlying C code, which
    leads to a major speed increase for large data sets.  Essentially,
    for the lower leaves all our time used to be eaten up by
    bookkeeping, and this was replaced by a different approach.
    The primary routine also uses .Call so as to be more memory
    efficient.

    \item The other major change was an error for asymmetric loss
    matrices, prompted by a user query.  
    With L=loss asymmetric, the altered priors were computed wrong -- they
were using L' instead of L.  Upshot -- the tree would not not necessarily
choose optimal splits for the given loss matrix.  Once chosen, splits were 
evaluated correctly.  The printed "improvement" values are of course the
wrong onese as well.  It's interesting that for my little test case, with
L quite asymmetric, the early splits in the tree are unchanged -- a good split
still looks good.
    \item Add the return.all argument to xpred.rpart.

    \item Added a set of formal tests, i.e., cases with known answers to
    which we can compare.
    \item Add a usercode vignette, explaining how to add user defined
    splitting functions.
    \item The class method now also returns the node probability.
    \item Add the stagec data set, used in some tests.
    \item The plot.rpart routine needs to store a value away that will
    be visible to the rpartco routine at a later time.  This is now
    stored in the options() list as rpartplot.
}}

\section{Changes in version 3.1-55, 2012 Oct 02 }{
  \itemize{
    \item Force use of registered symbols in R >= 2.16.0
    \item Update Polish translations.
    \item Work on message formats.
}}

\section{Changes in version 3.1-54, 2012 Jun 27 }{
  \itemize{
    \item Add Polish translations
}}

\section{Changes in version 3.1-53, 2012 Jun 01}{
  \itemize{
    \item rpart, rpart.matrix: allow backticks in formulae.
    \item tests/backtick.R: regession test
}}

\section{Changes in version 3.1-52, 2012 Mar 04}{
  \itemize{
    \item src/xval.c: ensure unused code is not compiled in.
}}

\section{Changes in version 3.1-51, 2012 Jan 11}{
  \itemize{
    \item  Change description of 'margin' in ?plot.rpart as suggested by
    Bill Venables.
}}    



