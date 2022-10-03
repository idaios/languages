#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "binary_tree.h"

#define SPACES 10

/**
 * Frees the resources associated to a binary tree with a recoursive depth-first visit
 *
 * @param node Root of a binary tree
 */
static void btree_free_df_r(pbtree_node node) {
    if (!node) return;

    btree_free_df_r(node->left);
    btree_free_df_r(node->right);

    free(node->name);
    memset(node, 0, sizeof(btree_node));
    free(node);
}

/**
 * Return the depth of a binary tree
 *
 * @param node Root of a binary tree
 * @returns Depth of a binary tree
 */
static size_t tree_get_depth(const btree_node * tree) {
    if (!tree) return 0;

    size_t left_depth = tree_get_depth(tree->left);
    size_t right_depth = tree_get_depth(tree->right);

    return 1 + ((left_depth > right_depth) ? left_depth : right_depth);
}

/**
 * Print recursively the binary tree to the console
 * @param tree Node that represent the oot of the tree at the current level
 * @param parent Parent of the current node
 * @param depth Level of depth for the current node.
 * @param right Indicates if we are in a right branch
 * @param levels Helper array
 */
static void
btree_draw_horizontal_r(const btree_node * tree, const btree_node * parent, const int depth, const int right,
                        char *levels) {
    if (tree == NULL) return;

    // Let's go with the right node first
    btree_draw_horizontal_r(tree->right, tree, depth + 1, 1, levels);

    if (depth > 0) {
        for (int i = 1; i < depth; i++) {
            // starting from level 1,
            printf(levels[i] ? "|" : " ");
            for (int j = 1; j < SPACES; j++)
                printf(" ");
        }
        // if we are printing a node that is not the root of the three,
        // we need to also print the edges with the associated value
        printf("+--%.2f-->", right ? parent->rvalue : parent->lvalue);
    }
    // Print the current node name
    printf("[%s]\n", tree->name);

    // levels keeps track of the direction of the branch that we have traversed at a certain depth
    // levels[i] is true if we have taken the left branch at level i to get to the current node
    levels[depth] = !levels[depth];
    // Let's go to left node
    btree_draw_horizontal_r(tree->left, tree, depth + 1, 0, levels);
}

void binary2trip( const btree_node *btree, ttree_node **ttree, ttree_node **rootpointer, char* deletedRootName, int restarti){
  static int i = 0;
  if(restarti) i = 0;
  i++;
  if(btree == NULL) return;

  if(strcmp(deletedRootName, btree->name) == 0){ return; }
 
  (*ttree) = calloc(1, sizeof(ttree_node));
  (*ttree)->name = calloc(strlen(btree->name)+1, sizeof(char));
  strcpy((*ttree)->name, btree->name);
  
  if(rootpointer == NULL)
    {
      fprintf(stderr, "ROOT IS NULL\n");
    }

  if(i == 1){
    *rootpointer = *ttree;
    assert(*rootpointer != NULL);

    if( ((btree)->left)->left != NULL && ((btree)->left)->right != NULL){
      //(btree->left)->delete = 1;
      strcpy(deletedRootName, (btree->left)->name);
      binary2trip((btree->left)->left, &((*ttree)->left), rootpointer, deletedRootName, 0);
      binary2trip((btree->left)->right, &((*ttree)->middle), rootpointer, deletedRootName, 0);
      
    }
    else if( ((btree)->right)->left != NULL && ((btree)->right)->right != NULL){
      //(btree->right)->delete = 1;
      strcpy(deletedRootName, (btree->right)->name);
      binary2trip( (btree->right)->left, &((*ttree)->right), rootpointer, deletedRootName, 0);
      binary2trip( (btree->right)->right, &((*ttree)->middle), rootpointer, deletedRootName, 0);
    }
    else{
      assert(0);
    }
  }
  else{
    (*ttree)->middle = NULL;
  }
    
  if( i != 1)
    {
      binary2trip(btree->left, &((*ttree)->left), rootpointer, deletedRootName, 0);
      binary2trip(btree->right, &((*ttree)->right), rootpointer, deletedRootName, 0);
    }
}


static void btree_draw_newick_b(const btree_node * tree, const btree_node * parent, const int depth, const int right,
			      char *levels, char* newickFormat) {
  if (tree == NULL) return;
  
  if(tree->left == NULL && tree->right == NULL){
    char *tmp = calloc(256, sizeof(char));
    sprintf(tmp, "%s:1", tree->name);
    strcat(newickFormat, tmp);
    free(tmp);
  }else{
    strcat(newickFormat, "(");
    btree_draw_newick_b(tree->left, tree, depth + 1, 1, levels, newickFormat);
    strcat(newickFormat, ",");
    btree_draw_newick_b(tree->right, tree, depth + 1, 1, levels, newickFormat);
    strcat(newickFormat, ")");
  }
}
  

static void ttree_draw_newick_b(const ttree_node * tree, char* newickFormat) {
  if (tree == NULL) return;
  
  if(tree->left == NULL && tree->right == NULL && tree->middle == NULL){
    char *tmp = calloc(256, sizeof(char));
    sprintf(tmp, "%s:1", tree->name);
    strcat(newickFormat, tmp);
    free(tmp);
  }else{
    strcat(newickFormat, "(");
    ttree_draw_newick_b(tree->left, newickFormat);
    strcat(newickFormat, ",");
    if(tree->middle != NULL)
      {
	ttree_draw_newick_b(tree->middle, newickFormat);
	strcat(newickFormat, ",");
      }
    ttree_draw_newick_b(tree->right, newickFormat);
    strcat(newickFormat, "):1");
  }
}

/**
 * Print a binary tree horizontally on the console
 *
 * @param tree Root of a binary tree
 */
void btree_draw_horizontal(const btree_node *tree) {
  if (!tree) return;
    // Let's use a support array in order to print the edges for the tree
    char *levels = calloc(tree_get_depth(tree), sizeof(char));
    btree_draw_horizontal_r(tree, NULL, 0, 0, levels);
    free(levels);
}

/**
 * Print a binary tree horizontally on the console
 *
 * @param tree Root of a binary tree
 */
void btree_draw_newick(const btree_node *tree, char* newicktree) {
  strcpy(newicktree, "");
  newicktree[1] = '\0';
  if (!tree) return;
  // Let's use a support array in order to print the edges for the tree
  char *levels = calloc(tree_get_depth(tree), sizeof(char));
  btree_draw_newick_b(tree, NULL, 0, 0, levels, newicktree);
  strcat(newicktree, ";");
  free(levels);
}



void ttree_draw_newick(const btree_node *btree, char* newicktree) {
  ttree_node *tree = NULL, *rootpointer = NULL;
  char *deletedRootName = calloc(256, sizeof(char));
  binary2trip(btree, &tree, &rootpointer, deletedRootName, 1);

  free(deletedRootName);
  tree = rootpointer;
  assert(rootpointer != NULL);
  
  if (!tree){
    fprintf(stderr, "Tree is null\n\n");
    return;
  }

    // Let's use a support array in order to print the edges for the tree
  ttree_draw_newick_b(tree, newicktree);
  newicktree[ strlen(newicktree) - 2 ] = ';';
  newicktree[ strlen(newicktree) - 1 ] = '\0';
}
 
/**
 * Frees the resources associated to a binary tree
 *
 * @param node Root of a binary tree
 */
 void btree_free(pbtree_node tree) {
   btree_free_df_r(tree);
 }
