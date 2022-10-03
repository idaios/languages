#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "fasta_parser.h"

struct _fasta_sequences {
    GSList *sequences; ///< List of all sequences read from a FASTA file
    size_t sequences_count; ///< Number of parsed sequences
};


/**
 * Substitutes the first occurrence of the characters [\n \r] with the \0 char
 *
 * @param str Source buffer
 */
static void line_replace_trailing(char *str) {
    char *pos;
    if ((pos = strchr(str, '\n')) != NULL)
        *pos = '\0';
    if ((pos = strchr(str, '\r')) != NULL)
        *pos = '\0';
}

/**
 * Parses a single FASTA sequence inside a file
 *
 * @param file_handle Handle of the opened file
 * @param line_buffer Buffer used for read the file (to avoid more allocations)
 * @param lineb_length Length of the buffer
 * @return Parsed FASTA sequence
 */
static char *fasta_parse_sequence(FILE *file_handle, char **line_buffer, size_t *lineb_length, size_t *sequence_length) {
    GSList *lines_list = NULL;
    ssize_t readed_chars = 0;
    *sequence_length = 0;
    char *buffer;

    //while ((readed_chars = getline(line_buffer, lineb_length, file_handle)) != -1) {
        // The buffer returner from getline contains the newline escape characters. We first need to remove them
    readed_chars = getline(line_buffer, lineb_length, file_handle);
    line_replace_trailing(*line_buffer);
    // We have encountered an empty line. The sequence is finished
    //if (strlen(*line_buffer) <= 0) break;
    // We need to copy the buffer before the next getline
    char *localb = calloc(readed_chars + 1, sizeof(char));
    strcpy(localb, *line_buffer);
    lines_list = g_slist_append(lines_list, localb);
    (*sequence_length) += readed_chars-1;
    
    // We must have read least one line of the sequence, otherwise the file is not valid
    assert(*sequence_length > 0);

    // We now only have to merge the read lines
    buffer = calloc(*sequence_length + 1, sizeof(char));
    for (GSList *item = lines_list; item != NULL; item = g_slist_next(item)) {
        strcat(buffer, item->data);
        free(item->data);
    }
    g_slist_free(lines_list);
    return buffer;
}

pfasta fasta_parse_file(const char *file_name, size_t *sequence_length, char*** fastaNames) {
    GSList *seq_list = NULL;
    char *curr_line = NULL;
    size_t line_length = 0;
    int InitialSize = 10000, i = 0, j = 0;
    int InitialLength = 256;
    char f;
    
    (*fastaNames) = calloc(InitialSize, sizeof(char*));
    for( i = 0; i < InitialSize; ++i){
      (*fastaNames)[i] = calloc(InitialLength, sizeof(char));
      assert( (*fastaNames)[i]  != NULL);
    }
    assert( fastaNames != NULL);
    assert( *fastaNames != NULL);
     
    FILE *file_handle = fopen(file_name, "r");
    if (!file_handle) return NULL;
    pfasta result = malloc(sizeof(fasta));
    if (!result) return NULL;
    
    i = 0;
    while (getline(&curr_line, &line_length, file_handle) != -1) {
        if (curr_line[0] == '>') {
	  
	  j = 0;
	  f = curr_line[1];
	  while(j < InitialLength && f != 13 && f!= 32 && f!= 9 && f!=10){
	    (*fastaNames)[i][j++] = f;
	    f = curr_line[j+1];
	  }
	  (*fastaNames)[i][j] = '\0';
	  i++;
            // We are reading a sequence
	  char *seq_buffer = fasta_parse_sequence(file_handle, &curr_line, &line_length, sequence_length);
            seq_list = g_slist_append(seq_list, seq_buffer);
        }
    }
    
    free(curr_line);
    fclose(file_handle);

    result->sequences = seq_list;
    result->sequences_count = g_slist_length(seq_list);
    return result;
}


void fasta_free(pfasta fasta_seq) {
    if (!fasta_seq) return;

    if (fasta_seq->sequences) {
        for (GSList *item = fasta_seq->sequences; item != NULL; item = g_slist_next(item)) {
            free(item->data);
        }
        g_slist_free(fasta_seq->sequences);
    }
    memset(fasta_seq, 0, sizeof(fasta));
    free(fasta_seq);
}

size_t fasta_sequences_count(const fasta *fasta_seq) {
    if (!fasta_seq) return 0;
    return fasta_seq->sequences_count;
}

GSList *fasta_sequences_list(const fasta *fasta_seq) {
    if (!fasta_seq) return 0;
    return fasta_seq->sequences;

}


