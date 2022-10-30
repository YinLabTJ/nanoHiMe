#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "training/mcg/kmer2binary.h"
#include "training/mcg/viterbi_algorithm.h"
#include "htslib/faidx.h"
#include <zlib.h>

int r_region_to_e_region(int E_loc[], int event_num, int r_start, int r_end, int * e_start, int * e_end, int start_site);

int model_calling(char model_file[], char chr[], char read_id[], char cref_seq[], int Eloc[], float Emean[], float Estd[], int event_num, int start_site, int rev_chain, char rev_mark, float Eduration[]);

int main(int argc, char *argv[])
{

	if(argc != 4) // ./viterbi model file, ref sequence and event_file
	{
		printf("Usage: re_align model_file ref_file event_file\n");
		return -1;
	}
	
	faidx_t *fai = fai_load(argv[2]);
	int fetch_len=0;

	int rev_chain; char rev_mark='+';
	int start_site,end_site;
	int event_num=0;

	float Emean[100000],Estd[100000]; int Eloc[100000]; float Eduration[100000];
	int e_total=0;

	float Emean_tmp,Estd_tmp,Eduration_tmp; int Eloc_tmp; int last_loc=0;
	char read_id[50]; char last_read_id[50]="None"; 
	char chr[30]; char ref_kmer[7],model_kmer[7];
	int pos_c=0; int neg_c=0;

	char * cref_seq;


	gzFile fi;
	fi = gzopen(argv[3],"rb");
	char string[1000];
	while(gzgets(fi, string, 1000) != NULL){
		sscanf(string,"%s %d %s %s %*s %*s %f %f %f %s %*s %*s %*s",chr,&Eloc_tmp,ref_kmer,read_id,&Emean_tmp,&Estd_tmp,&Eduration_tmp,model_kmer);

		if(strcmp(chr, "contig")!=0){
											
			if((strcmp(read_id, last_read_id) == 0 && Eloc_tmp-last_loc<=20)|| strcmp(last_read_id, "None")==0){
				Eloc[e_total]=Eloc_tmp;
				last_loc=Eloc_tmp;
				Emean[e_total]=Emean_tmp;
				Estd[e_total]=Estd_tmp;
				Eduration[e_total]=Eduration_tmp;
			}else{
				start_site = Eloc[0]; end_site = Eloc[e_total-1]+5;
				cref_seq = faidx_fetch_seq(fai, chr , start_site , end_site , &fetch_len);					
				if(pos_c>neg_c){
					rev_chain=0; rev_mark='+';
				}else{
					rev_chain=1; rev_mark='-';
				}
				model_calling(argv[1], chr, last_read_id, cref_seq, Eloc, Emean, Estd, e_total, start_site, rev_chain, rev_mark, Eduration);
			
				e_total=0; pos_c=0; neg_c=0;

				Eloc[0]=Eloc_tmp;
				last_loc=Eloc_tmp;
				Emean[0]=Emean_tmp;
				Estd[0]=Estd_tmp;
				Eduration[0]=Eduration_tmp;
			}
			e_total++;
			strcpy(last_read_id, read_id);
			if(model_kmer == "NNNNNN"){
			}else{
				if(strcmp(ref_kmer, model_kmer) == 0){
					pos_c++;
				}else{
					neg_c++;
				}
			}
		}
	}

	start_site = Eloc[0]; end_site = Eloc[e_total-1]+5;
	cref_seq = faidx_fetch_seq(fai, chr , start_site , end_site , &fetch_len);
	if(pos_c>neg_c){
		rev_chain=0; rev_mark='+';
	}else{
		rev_chain=1; rev_mark='-';
	}
	model_calling(argv[1], chr, read_id, cref_seq, Eloc, Emean, Estd, e_total, start_site, rev_chain, rev_mark, Eduration);

	return 0;
}


int model_calling(char model_file[], char chr[], char read_id[], char cref_seq[], int Eloc[], float Emean[], float Estd[], int event_num, int start_site, int rev_chain, char rev_mark, float Eduration[]){
	float score;
	score = viterbi_score(model_file,cref_seq,Eloc,Emean,Estd,chr,start_site,event_num,rev_chain,0,read_id,Eduration);

	return 0;
}

int r_region_to_e_region(int E_loc[], int event_num, int r_start, int r_end, int * e_start, int * e_end, int start_site){
	
	*e_start=0; *e_end=0;
	for(long i=0; i<event_num; i++){
		if(E_loc[i]>=r_start+start_site && E_loc[i]<=r_end+start_site){
			if(*e_start==0)
				*e_start=i;
			if(*e_end==0)
				*e_end=i;
			if(*e_start>i)
				*e_start=i;
			if(*e_end<i)
				*e_end=i;
		}
	}
	return 0;
}
