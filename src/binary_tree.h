#ifndef _BINARY_TREE_H
#define _BINARY_TREE_H

#include <assert.h>

typedef struct _btree_node btree_node, *pbtree_node;
typedef struct _ttree_node ttree_node, *pttree_node;

typedef struct _btree_node {
  int delete;
  char *name; ///< Value of the node
  pbtree_node left; ///< Left node
  double lvalue; ///< Value of the edge to the left node
  pbtree_node right; ///< Right node
  double rvalue; ///< Value of the edge to the right node
} btree_node, *pbtree_node;


typedef struct _ttree_node {
  char *name; ///< Value of the node
  pttree_node left; ///< Left node
  double lvalue; ///< Value of the edge to the left node
  pttree_node right; ///< Right node
  double rvalue; ///< Value of the edge to the right node
  pttree_node middle;
  double mvalue;
} ttree_node, *pttree_node;

void btree_draw_horizontal(const btree_node *tree);
void btree_draw_newick(const btree_node *tree, char* newick);
void ttree_draw_newick(const btree_node *tree, char* newick);
void btree_free(pbtree_node tree);

#endif // !_BINARY_TREE_H
