/*
** SCCS @(#)rpartS.h	1.2 02/08/98
** some macros that are needed for those routines that call "S_alloc"
**  they aren't included "straigtht" into the routines because they differ
**  between Splus and S4 (Bell Labs 4, which will be Splus 5)
*/

/* S4 version */
/* #include "S.h"
** #define ALLOC(a,b) S_alloc(a,b,S_evaluator)
*/

/* Splus version */
#define ALLOC(a,b) S_alloc(a,b)

extern char *S_alloc();  
