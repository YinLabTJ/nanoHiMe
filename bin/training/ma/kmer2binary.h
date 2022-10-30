#include <string.h>
#include <math.h>
int kmer_to_num(char * kmer);

int kmer_to_num_cg(char * kmer);

void reverse_complement(char nput_seq[], char * output_seq);

void methylation_transition(char input_seq[], char * output_seq);

void methylation_transition_neg(char input_seq[], char * output_seq);

float includeA(char input_seq[]);
