#include <string.h>
#include <math.h>
int kmer_to_num(char * kmer);

void reverse_complement(char nput_seq[], char * output_seq);

void methylation_transition(char input_seq[], char * output_seq);

void methylationA_transition(char input_seq[], char * output_seq);

void methylation_transition_neg(char input_seq[], char * output_seq);

void methylationA_transition_neg(char input_seq[], char * output_seq);
