#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "ma_region/kmer2binary.h"
#include "ma_region/viterbi_algorithm.h"
#include "htslib/faidx.h"
#include <zlib.h>

int r_region_to_e_region(int E_loc[], int event_num, int r_start, int r_end, int * e_start, int * e_end, int start_site, int * r_start_t, int * r_end_t);

int ma_region_calling(char model_dir[], char chr[], char read_id[], char cref_seq[], int Eloc[], float Emean[], float Estd[], int event_num, int start_site, int rev_chain, char rev_mark, int win_width, int step_width);

int main(int argc, char *argv[])
{

	if(argc != 6) 
	{
		printf("Usage: methylate_ma model_dir ref_file event_file win_width step_width\n");
		return -1;
	}
	
	faidx_t *fai = fai_load(argv[2]);
	int fetch_len=0;

	int win_width=atoi(argv[4]);
	int step_width=atoi(argv[5]);
	int rev_chain; char rev_mark='+';
	int start_site,end_site;
	int event_num=0;

	

	float Emean[100000],Estd[100000]; int Eloc[100000];
	int e_total=0;
	int event_num_t;

	float Emean_tmp,Estd_tmp; int Eloc_tmp;
	char read_id[50]; char last_read_id[50]="None"; 
	char chr[30],last_chr[30]; char ref_kmer[7],model_kmer[7];
	int pos_c=0; int neg_c=0;

	char * cref_seq;
	
	gzFile fi;
	fi = gzopen(argv[3],"rb");
	char string[1000];
	while(gzgets(fi, string, 1000) != NULL){
		sscanf(string,"%s %d %s %s %*s %*s %f %f %*s %s %*s %*s %*s",chr,&Eloc_tmp,ref_kmer,read_id,&Emean_tmp,&Estd_tmp,model_kmer);
		
		if(strcmp(chr, "contig")!=0){
			if(strcmp(read_id, last_read_id) == 0 || strcmp(last_read_id, "None")==0){
				Eloc[e_total]=Eloc_tmp;
				Emean[e_total]=Emean_tmp;
				Estd[e_total]=Estd_tmp;
			}else{
				start_site = Eloc[0]; end_site = Eloc[e_total-1]+5;
				cref_seq = faidx_fetch_seq(fai, last_chr , start_site , end_site , &fetch_len);					
				if(pos_c>neg_c){
					rev_chain=0; rev_mark='+';
				}else{
					rev_chain=1; rev_mark='-';
				}
				ma_region_calling(argv[1], last_chr, last_read_id, cref_seq, Eloc, Emean, Estd, e_total, start_site, rev_chain, rev_mark, win_width, step_width);
			
				e_total=0; pos_c=0; neg_c=0;

				Eloc[0]=Eloc_tmp;
				Emean[0]=Emean_tmp;
				Estd[0]=Estd_tmp;
			}
			e_total++;
			strcpy(last_read_id, read_id);
			strcpy(last_chr, chr);
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
	ma_region_calling(argv[1], chr, read_id, cref_seq, Eloc, Emean, Estd, e_total, start_site, rev_chain, rev_mark, win_width, step_width);

	return 0;
}


int ma_region_calling(char model_dir[],char chr[], char read_id[], char cref_seq[], int Eloc[], float Emean[], float Estd[], int event_num, int start_site, int rev_chain, char rev_mark, int win_width, int step_width){
	int len=strlen(cref_seq);
	
	char model_p[200],model_a[200],model_cg[200],model_acg[200];
	strcpy(model_p,model_dir); strcat(model_p,"/PCR/PCR.param.txt");
	strcpy(model_a,model_dir); strcat(model_a,"/6mA/6mA.param.txt");
	strcpy(model_cg,model_dir); strcat(model_cg,"/mCG/mCG.param.txt");
	strcpy(model_acg,model_dir); strcat(model_acg,"/6mA_mCG/6mA_mCG.param.txt");
	

	float unmethylscore,cgmethylscore,amethylscore1,amethylscore2,acgmethylscore1,acgmethylscore2;
	
	int r_start=0; int r_end=0; int e_start; int e_end; int r_start_t=0; int r_end_t=0;
	char call_region[1500]; 
	for(int i=0; i<len-win_width; i+=step_width){
			r_start=i; r_end=i+win_width-1; 
			if(r_end>len)
				r_end=len;
			r_region_to_e_region(Eloc, event_num, r_start,r_end-5, &e_start, &e_end, start_site, &r_start_t, &r_end_t);
			if(e_start==0 && e_end==0){
				strcpy(call_region,"NNNNNNNNNN");
			}else{
				strncpy(call_region,&cref_seq[r_start_t],r_end_t-r_start_t+1);
				call_region[r_end_t-r_start_t+1]='\0';
			}

			int e_loc[e_end-e_start+1];
			float e_mean[e_end-e_start+1],e_std[e_end-e_start+1];
			for(int k=0;k<=e_end-e_start;k++)
			{
				e_loc[k]=Eloc[e_start+k];
				e_mean[k]=Emean[e_start+k];
				e_std[k]=Estd[e_start+k];
			}
			
			unmethylscore=viterbi_score(model_p,call_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
			cgmethylscore=viterbi_score(model_cg,call_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
			amethylscore1=viterbi_score(model_a,call_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
			amethylscore2=viterbi_score(model_a,call_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,1);
			acgmethylscore1=viterbi_score(model_acg,call_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
			acgmethylscore2=viterbi_score(model_acg,call_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,1);
			float amethylscore=amethylscore1;
			if(amethylscore<amethylscore2)
				amethylscore=amethylscore2;
			float acgmethylscore=acgmethylscore1;
			if(acgmethylscore<acgmethylscore2)
				acgmethylscore=acgmethylscore2;
			printf("%s\t%c\t%d\t%d\t%s\t%.2f\t%s\n",chr,rev_mark,start_site+r_start,start_site+r_end,read_id,acgmethylscore-cgmethylscore+amethylscore-unmethylscore,call_region);
			//printf("%s\t%c\t%d\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",chr,rev_mark,start_site+r_start,start_site+r_end,read_id,acgmethylscore-cgmethylscore+amethylscore-unmethylscore,acgmethylscore-amethylscore+cgmethylscore-unmethylscore,acgmethylscore1,acgmethylscore2,amethylscore1,amethylscore2,cgmethylscore,unmethylscore,call_region);
	}

	return 0;
}

int r_region_to_e_region(int E_loc[], int event_num, int r_start, int r_end, int * e_start, int * e_end, int start_site, int * r_start_t, int * r_end_t){
	
	*e_start=0; *e_end=0;
	for(long i=0; i<event_num; i++){
		if(E_loc[i]>=r_start+start_site && E_loc[i]<=r_end+start_site){
			if(*e_start==0){
				*e_start=i; 
				*r_start_t=E_loc[i]-start_site;
			}
			if(*e_end==0){
				*e_end=i;
				*r_end_t=E_loc[i]+5-start_site;
			}
			if(*e_start>i){
				*e_start=i;
				*r_start_t=E_loc[i]-start_site;
			}
			if(*e_end<i){
				*e_end=i;
				*r_end_t=E_loc[i]+5-start_site;
			}
		}
	}
	return 0;
}
