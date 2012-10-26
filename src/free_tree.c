/*
** free up all of the memory associated with a tree
*/
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

static void free_split(struct split *spl) 
{
    if (spl != 0) {
	free_split(spl->nextsplit);
	Free(spl);
    }
}
void free_tree(struct node *node,  int freenode)
{
    if (node->rightson != 0) free_tree(node->rightson, 1);
    if (node->leftson  != 0) free_tree(node->leftson,  1);

    free_split(node->surrogate);
    free_split(node->primary);
    if (freenode == 1) Free(node);
}


	    
