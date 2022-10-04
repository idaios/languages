#include<stdio.h>
#include<stdlib.h>
#include <levenshtein.h>
#include <assert.h>
#include "clusters_matrix.h"
#include "fasta_parser.h"
#include "cmdline.h"
#include <time.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define NULLCHAR '?'
#include <libpll/pll_tree.h>
#include <assert.h>
#include <stdarg.h>

/* set to 1 for printing splits */
#define PRINT_SPLITS 0

/* static functions */
static void fatal (const char * format, ...);

int rf(char string1[], char string2[])
{
  unsigned int rf_dist;

  /* tree properties */
  pll_utree_t * tree1 = NULL,
              * tree2 = NULL;
  pll_rtree_t * tree1rooted = NULL,* tree2rooted=NULL;
  unsigned int tip_count;

  /* if (argc != 3) */
  /*   fatal (" syntax: %s [newick] [newick]", argv[0]); */

  
  /* parse the input trees */
  tree1 = pll_utree_parse_newick_string(string1);
  tree2 = pll_utree_parse_newick_string(string2);
  tip_count = tree1->tip_count;
  
  if (tip_count != tree2->tip_count)
    {
      fprintf(stderr, "Trees have different number of tips %d and %d\n", tip_count, tree2->tip_count);
      fatal("Trees have different number of tips!");
    }

  if (!pllmod_utree_consistency_set(tree1, tree2))
    fatal("Cannot set trees consistent!");

  if (!pllmod_utree_consistency_check(tree1, tree2))
    fatal("Tip node IDs are not consistent!");



  /* uncomment lines below for displaying the trees in ASCII format */
  //pll_utree_show_ascii(tree1, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);
  //pll_utree_show_ascii(tree2, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

  /* next, we compute the RF distance in 2 different ways: */

  /* 1. creating the split sets manually */
  unsigned int n_splits = tip_count - 3;
  pll_split_t * splits1 = pllmod_utree_split_create(tree1->nodes[tip_count],
                                                    tip_count,
                                                    NULL);

#if(PRINT_SPLITS)
  {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits1[i], tip_count);
      printf("\n");
    }
    printf("\n");
  }
#endif

  /* now we compute the splits, but also the nodes corresponding to each split */
  pll_unode_t ** splits_to_node = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits2 = pllmod_utree_split_create(tree2->nodes[tip_count],
                                                    tip_count,
                                                    splits_to_node);

#if(PRINT_SPLITS)
  {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits2[i], tip_count);
      printf(" node: Pmatrix:%d Nodes:%d<->%d\n",
                                splits_to_node[i]->pmatrix_index,
                                splits_to_node[i]->node_index,
                                splits_to_node[i]->back->node_index);
    }
  }
#endif

  rf_dist = pllmod_utree_split_rf_distance(splits1, splits2, tip_count);
  printf("RF [manual]\n");
  printf("distance = %d\n", rf_dist);
  printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  pllmod_utree_split_destroy(splits1);
  pllmod_utree_split_destroy(splits2);
  free(splits_to_node);

  /* 2. directly from the tree structures */

  rf_dist = pllmod_utree_rf_distance(tree1->nodes[tip_count],
                                     tree2->nodes[tip_count],
                                     tip_count);

  printf("RF [auto]\n");
  printf("distance = %d\n", rf_dist);
  printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  /* clean */
  pll_utree_destroy (tree1, NULL);
  pll_utree_destroy (tree2, NULL);

  return (0);
}



double rf2(char string1[], pll_utree_t *target_tree)
{
  unsigned int rf_dist;

  /* tree properties */
  pll_utree_t * tree1 = NULL,
              * tree2 = NULL;
  pll_rtree_t * tree1rooted = NULL,* tree2rooted=NULL;
  unsigned int tip_count;
  
  /* parse the input trees */
  //  fprintf(stderr, "%s\n", string1);
  tree1 = pll_utree_parse_newick_string(string1);
  tree2 = target_tree; //pll_utree_parse_newick_string(string2);
  
  tip_count = tree2->tip_count;
    
  if (tip_count != tree1->tip_count)
    {
      fprintf(stderr, "Trees have different number of tips %d and %d\n", tip_count, tree2->tip_count);
      fatal("Trees have different number of tips!");
    }

  if (!pllmod_utree_consistency_set(tree1, tree2))
    fatal("Cannot set trees consistent!");

  if (!pllmod_utree_consistency_check(tree1, tree2))
    fatal("Tip node IDs are not consistent!");

  //fprintf(stderr, "tip count: %d\n", tip_count); 

  /* uncomment lines below for displaying the trees in ASCII format */
  //pll_utree_show_ascii(tree1, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);
  //pll_utree_show_ascii(tree2, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

  /* next, we compute the RF distance in 2 different ways: */

  /* 1. creating the split sets manually */
  unsigned int n_splits = tip_count - 3;
  pll_split_t * splits1 = pllmod_utree_split_create(tree1->nodes[tip_count],
                                                    tip_count,
                                                    NULL);

  /* now we compute the splits, but also the nodes corresponding to each split */
  pll_unode_t ** splits_to_node = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits2 = pllmod_utree_split_create(tree2->nodes[tip_count],tip_count,splits_to_node);
  
  rf_dist = pllmod_utree_split_rf_distance(splits1, splits2, tip_count);
  
  /* printf("RF [manual]\n"); */
  /* printf("distance = %d\n", rf_dist); */
  /* printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3))); */

  /* pllmod_utree_split_destroy(splits1); */
  /* pllmod_utree_split_destroy(splits2); */
  /* free(splits_to_node); */

  /* /\* 2. directly from the tree structures *\/ */

  /* rf_dist = pllmod_utree_rf_distance(tree1->nodes[tip_count], */
  /*                                    tree2->nodes[tip_count], */
  /*                                    tip_count); */

  /* printf("RF [auto]\n"); */
  /* printf("distance = %d\n", rf_dist); */
  /* printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3))); */

  /* clean */
  //pll_utree_destroy (tree1, NULL);
  //pll_utree_destroy (tree2, NULL);

  //return (rf_dist/(2.0*(tip_count-3)));
  return(rf_dist);

  return (0);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}

/////////////////////




int hamming(char **seqs, int nseqs, size_t length, int*** distMatrix, int*** LengthMatrix){
  int d = 0;
  size_t i = 0;
  int j1, j2; 
  
  for(i =0; i < length; ++i)
    {
      for(j1 = 0; j1 < nseqs; ++j1){
	if(seqs[j1][i] != NULLCHAR)
	  LengthMatrix[i][j1][j1] = 1;
	for(j2 = j1+1; j2 <nseqs; ++j2){
	  if( seqs[j1][i] != NULLCHAR && seqs[j2][i] != NULLCHAR){
	    LengthMatrix[i][j1][j2] = 1;
	    LengthMatrix[i][j2][j1] = 1;
	    if( seqs[j1][i] != seqs[j2][i] ){
	      distMatrix[i][j1][j2] = distMatrix[i][j2][j1] = 1;
	    }
	  }
	}
      }
    }
  return 1;
}

int getNJDistance(int* njdistanceMatrix, int* validComparisons, int k, int* indexes, int*** matrices, int*** validComparisonsMatrices, int sz){
  int i = 0, ii = 0, iii = 0;

  for(i = 0; i < sz*sz; ++i){
    njdistanceMatrix[i] = 0;
    validComparisons[i] = 0;
  }
  
  for(i = 0; i <k; ++i){
    for(ii = 0; ii < sz; ++ii){
      for(iii = 0; iii < sz; ++iii){
	njdistanceMatrix[ii*sz + iii] += matrices[ indexes[i] ][ii][iii];
	validComparisons[ii*sz + iii] += validComparisonsMatrices[ indexes[i] ][ii][iii];
      }
    }
  }
  return 1;
}

void shuffleInts(int *array, int n)
{
  for(int i = 0; i < n; ++i)
    {
      int tmp = array[i];
      int randomIndex = rand() % n;
      array[i] = array[randomIndex];
      array[randomIndex] = tmp;
    }
}

void  binomial_random(int sz, double prob, int* inds, int* indssz)
{
  int cnt = 0;
  double rnd = 0.;
  for(int i = 0; i < sz; ++i)
    {
      rnd = ((double)rand())/RAND_MAX;
      if( rnd < prob)
	{
	  inds[cnt] = i;
	  cnt++;
	}
    }
  *indssz = cnt;
}

double subsetRF(int*** matrices, int* njdistanceMatrix, int k, int sz, int seqLength, int*** validComparisonsMatrices, int* validComparisons, char* newicktree, pll_utree_t *target_tree, char **fastaNames, int *currentChosenIndexes, int *indexesToChooseFrom, int *k1, int *k2, double prob, int *indexes, int *nextToTryIndexes, int*njrescaledMat, int* keepinds, int* insertinds, int *tmpFlags, double rfPrevious){
  double rf = 0;
  int i = 0, keepindssz = 0, insertindssz = 0;
  double probdelete = 0., probinsert = 0.;

  static int statici = 0;
  statici++;
  
  
  if(statici == 1){
    assert(k > 0 && k < seqLength);
    int* v = calloc(seqLength, sizeof(int));

    for(int i = 0; i<seqLength; ++i){
      v[i] = i;
    }
    shuffleInts(v, seqLength);
    
    for(int i = 0; i < k; ++i)
      {
	currentChosenIndexes[i] = indexes[i] = v[i];
	//fprintf(stderr, "%d,", v[i]);
      }
    *k1 = k;
    *k2 = seqLength - k;
    for(int i = 0; i < (*k2); ++i)
      {
	assert( i + (*k1) < seqLength);
	indexesToChooseFrom[i] = v[i+(*k1)];
      }
    free(v);
  }
  else
    {
      assert(prob < 1 && prob > 0);
      probdelete = 1. - prob;
      probinsert = probdelete * ( (double)(*k1)/(*k2) );
      binomial_random( (*k1), prob, keepinds, &keepindssz);
      binomial_random( (*k2), probinsert, insertinds, &insertindssz);

      //fprintf(stderr, "prob: %f, probinsert: %f, k1: %d, keep: %d and insert: %d\n", prob, probinsert, *k1, keepindssz, insertindssz);
      for(int i = 0; i < keepindssz; ++i)
	{
	  nextToTryIndexes[i] = currentChosenIndexes[ keepinds[i] ];
	  //fprintf(stderr, "next at position %d is %d\n", i, nextToTryIndexes[i]);
	}
      
      for(int i = 0; i < insertindssz; ++i)
	{
	  nextToTryIndexes[i+keepindssz] = indexesToChooseFrom[ insertinds[i] ];
	  //fprintf(stderr, "next at position %d is %d\n", i, nextToTryIndexes[i]);
	}
      k = keepindssz + insertindssz;
    }
  
  getNJDistance(njdistanceMatrix, validComparisons, k, nextToTryIndexes, matrices, validComparisonsMatrices, sz);

  for(int i = 0; i < sz*sz; ++i)
    {
      njrescaledMat[i] = (int)((njdistanceMatrix[i]*10000.0)/validComparisons[i]); // just a small trick to normalize with the number of comparisons, because downstread function need an integer distance matrix
      //printf("%d,", njrescaledMat[i]);
    }
  //printf("\n");
  
  pclusters_matrix matrix = clusters_create_matrix(njrescaledMat, sz, 0, fastaNames);
  
  if (!matrix) {
    perror("Unable to create clustering instance");
    assert(0);
  }

  while (clusters_get_count(matrix) > 1) {
    clusters_increase_clustering(matrix);
  }
  
  ttree_draw_newick(clusters_get_tree(matrix), newicktree);
  rf = rf2(newicktree, target_tree);


  if(statici == 1){
    return rf;
  }

  if(rf <=  rfPrevious)
    {

      if(rf < rfPrevious){
	fprintf(stdout, "step: %d, rf: %f, k1: %d\n", statici, rf, *k1);
      }
	
      
      (*k1) = (*k2) = 0;
      for(int i = 0; i < seqLength; ++i)
	{
	  tmpFlags[i] = 0;
	}
      for(int i=0; i < k; ++i)
	{
	  currentChosenIndexes[i] = nextToTryIndexes[i];
	  tmpFlags[ currentChosenIndexes[i] ] = 1;
	  //fprintf(stderr, "Index %d became 1\n", currentChosenIndexes[i]);
	  (*k1)++;
	}
      int j = 0;
      for(int i = 0; i < seqLength; ++i)
	{
	  if(tmpFlags[i] == 0){
	    indexesToChooseFrom[j++] = i;
	    (*k2)++;
	  }
	}
      if( (*k1) + (*k2) != seqLength){
	//fprintf(stderr, "k1: %d, k2: %d, seqLength: %d\n", *k1, *k2, seqLength);
	assert( (*k1) + (*k2) == seqLength );
      }
      return rf;
    }
  else{
    return rfPrevious;
  }
  
  return -1;
}

int main(int argc, char **argv) {
  srand ( time(NULL) );
  size_t i=0, j=0;
  size_t sequence_length = 0;
  int target_is_rooted = 1;
  double keepProb = 0.99;
  static struct gengetopt_args_info args_info;
  assert(cmdline_parser(argc, argv, &args_info) == 0);
  int verbose = args_info.verbose_flag;
  if (args_info.inputs_num == 0) {
    fprintf(stderr, "FASTA input file not specified");
    return -1;
  }
  
  char *source_file = args_info.inputs[0];
  char *target_tree_file = args_info.inputs[1];
  char *rootedUnrooted = args_info.inputs[2];
  if( strcmp(rootedUnrooted, "rooted") == 0){ target_is_rooted = 1; }
  else{ target_is_rooted = 0; }
  char **fastaNames = NULL; 
  
  pfasta parsed = fasta_parse_file(source_file, &sequence_length, &fastaNames);
  
  if (!parsed || fasta_sequences_count(parsed) == 0) {
    fasta_free(parsed);
    perror("Unable to parse FASTA file");
    return -1;
  }
  
  pll_utree_t *target_tree = NULL;
  if( target_is_rooted ){
    pll_rtree_t *target_tree_rooted = pll_rtree_parse_newick(target_tree_file);
    target_tree = pll_rtree_unroot(target_tree_rooted);
  }
  else{
    target_tree = pll_utree_parse_newick(target_tree_file);
  }
  
  size_t count = fasta_sequences_count(parsed);
  
  int ***distances_matrix = calloc(sequence_length, sizeof(int**));
  int ***validcomparisons_matrix = calloc(sequence_length, sizeof(int**));
  for(i = 0; i < sequence_length; ++i){
    distances_matrix[i] = calloc(count, sizeof(int*));
    validcomparisons_matrix[i] = calloc(count, sizeof(int*));
    for(j = 0; j < count; ++j){
      distances_matrix[i][j] = calloc(count, sizeof(int));
      validcomparisons_matrix[i][j] = calloc(count, sizeof(int));
    }
  }
  
  char** seqs = calloc(count, sizeof(char*));
  for(j = 0; j < count; ++j)
    seqs[j] = calloc(sequence_length+1, sizeof(char));
  
  j = 0;
  for (GSList *item = fasta_sequences_list(parsed); item != NULL; item = g_slist_next(item)) {
    assert(j < count);
    strcpy(seqs[j], item->data);
    ++j;
  }
  fasta_free(parsed);
  
  hamming(seqs, count, sequence_length, distances_matrix, validcomparisons_matrix);
  
  int k = 200;
  int* njdistanceMatrix = calloc(count*count, sizeof(int));
  int* validComparisons = calloc(count*count, sizeof(int));
  
  
  //getNJDistance(njdistanceMatrix, k, indexes, distances_matrix, count);
  //////////////////////////////////////////////////////////////////////////
  //void subsetRF(int*** matrices, int* njdistanceMatrix, int k, int sz, int seqLength, int** validComparisons, char* newicktree, pll_utree_t *target_tree, char **fastaNames){
  
  int *njrescaledMat = calloc(count*count, sizeof(int));
  int *indexesToChooseFrom = calloc(sequence_length, sizeof(int));
  int *currentChosenIndexes = calloc(sequence_length, sizeof(int));
  int *tmpFlags = calloc(sequence_length, sizeof(int));
  int *nextToTryIndexes = calloc(sequence_length, sizeof(int));
  int *keepinds = calloc(sequence_length, sizeof(int));
  int *insertinds = calloc(sequence_length, sizeof(int));
  int *currentBestIndexes = calloc(sequence_length, sizeof(int));
  double rfdist = 0;
  int k1 = 0, k2 = 0;
  double rfPrevious = 9999999999;
  for( int reps = 0; reps < 100000; ++reps)
    {
      char *newicktree = calloc(10000, sizeof(char));
      rfPrevious = subsetRF(distances_matrix, njdistanceMatrix, k, count, sequence_length, validcomparisons_matrix, validComparisons, newicktree, target_tree, fastaNames,  currentChosenIndexes, indexesToChooseFrom, &k1, &k2, keepProb, currentBestIndexes, nextToTryIndexes, njrescaledMat, keepinds, insertinds, tmpFlags, rfPrevious);
      //fprintf(stdout, "%d, rf: %f\n", reps, rfPrevious);
      free(newicktree);
    }
  
  //////////////////////////////////////////////////////////////////////////
  
  /* pclusters_matrix matrix = clusters_create_matrix(njdistanceMatrix, count, verbose, fastaNames); */
  
  /* if (!matrix) { */
  /*     free(distances_matrix); */
  /*     perror("Unable to create clustering instance"); */
  /*     return -1; */
  /* } */
  /* while (clusters_get_count(matrix) > 1) { */
  /*     clusters_increase_clustering(matrix); */
  /* } */
  
  
  
  
    /* char *newicktree = calloc(10000, sizeof(char)); */
    /* ttree_draw_newick(clusters_get_tree(matrix), newicktree); */
    /* //    fprintf(stderr, "newick\n%s\n", newicktree); */
    /* rf2(newicktree, target_tree); */
    /* clusters_free_distances(matrix); */
    /* free(distances_matrix); */
}
