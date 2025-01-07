/*
 * free up all of the memory associated with a tree
 */
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

static void
free_split(pSplit spl)
{
    if (spl) {
	free_split(spl->nextsplit);
	R_Free(spl);
    }
}

/* use freenode if the tree was CALLOC-ed, from xval.c */
void
free_tree(pNode node, int freenode)
{
    if (node->rightson)
	free_tree(node->rightson, 1);
    if (node->leftson)
	free_tree(node->leftson, 1);

    free_split(node->surrogate);
    free_split(node->primary);
    if (freenode == 1)
	R_Free(node);
    else {
       /* don't point to things I just freed */
	node->primary = (pSplit) NULL;
	node->surrogate = (pSplit) NULL;
	node->rightson = (pNode) NULL;
	node->leftson = (pNode) NULL;
    }
}
