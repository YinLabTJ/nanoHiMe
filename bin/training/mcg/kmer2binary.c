#include <stdio.h>
#include <math.h>
#include "kmer2binary.h"

int kmer_to_num(char * kmer)
{
	int num=0;
	for(int i=5;i>=0;i--){
		switch(kmer[i]){
			case 'C': num += 1 * pow(4,5-i); break;
			case 'G': num += 2 * pow(4,5-i); break;
			case 'T': num += 3 * pow(4,5-i); break;
			case 'A': break;
		}	
	}
	return num;
}

int kmer_to_num_cg(char * kmer)
{
	int num=0;
	for(int i=5;i>=0;i--){
		switch(kmer[i]){
			case 'C': num += 1 * pow(5,5-i); break;
			case 'G': num += 2 * pow(5,5-i); break;
			case 'M': num += 3 * pow(5,5-i); break;
			case 'T': num += 4 * pow(5,5-i); break;
			case 'A': break;
		}
	}
}


void reverse_complement(char input_seq[], char * output_seq)
{
	int len = strlen(input_seq);
	for(int i=len;i>0;i--){
		switch(input_seq[i-1]){
			case 'C': output_seq[len-i]='G'; break;
			case 'G': output_seq[len-i]='C'; break;
			case 'T': output_seq[len-i]='A'; break;
			case 'A': output_seq[len-i]='T'; break;
			case 'W': output_seq[len-i]='M'; break; //means the G of mCG on other chain
		}
	}
	output_seq[len]='\0';
}

void methylation_transition(char input_seq[], char * output_seq)
{
	strcpy(output_seq,input_seq);
	int len = strlen(output_seq);
	for(int i=0;i<len-1;i++){
		if(output_seq[i]=='C' && output_seq[i+1]=='G')
			output_seq[i]='M';
	}
}

void methylation_transition_neg(char input_seq[], char * output_seq)
{
	strcpy(output_seq,input_seq);
	int len = strlen(output_seq);
	for(int i=0;i<len-1;i++){
		if(output_seq[i]=='C' && output_seq[i+1]=='G')
			output_seq[i+1]='W';
	}
}

float includeA(char input_seq[])
{
	int len = strlen(input_seq);
	float incluA=0;
	for(int i=0;i<len;i++){
		if(input_seq[i]=='A')
			incluA=1.5;
	}
	return incluA;
}

float includeCG(char input_seq[])
{
	int len = strlen(input_seq);
	float incluCG=0;
	for(int i=0;i<len-1;i++){
		if(input_seq[i]=='C' && input_seq[i+1]=='G')
			incluCG=1;
	}
	return incluCG;
}
