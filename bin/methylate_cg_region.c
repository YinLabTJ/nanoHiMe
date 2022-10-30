#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "mcg_region/kmer2binary_cg.h"
#include "mcg_region/viterbi_algorithm_cg.h"
#include "htslib/faidx.h"
#include <zlib.h>

int r_region_to_e_region(int E_loc[], int event_num, int r_start, int r_end, int * e_start, int * e_end, int start_site, int * r_start_t, int * r_end_t);

int cg_region_calling(char model_dir[], char chr[], char read_id[], char cref_seq[], int Eloc[], float Emean[], float Estd[], int event_num, int start_site, int rev_chain, char rev_mark);

int main(int argc, char *argv[])
{

	if(argc != 4) // ./main model file, ref sequence and event_file
	{
		printf("Usage: methylate_cg model_dir ref_file event_file\n");
		return -1;
	}
	
	faidx_t *fai = fai_load(argv[2]);
	int fetch_len=0;

	int rev_chain; char rev_mark='+';
	int start_site,end_site;
	int event_num=0;


	float Emean[100000],Estd[100000]; int Eloc[100000];
	int e_total=0;

	float Emean_tmp,Estd_tmp; int Eloc_tmp;
	char read_id[50]; char last_read_id[50]="None"; 
	char chr[30]; char ref_kmer[7],model_kmer[7];
	int pos_c=0; int neg_c=0;

	char * cref_seq;
	
	gzFile fi;
	fi = gzopen(argv[3],"rb");
	char string[1000];
	gzgets(fi, string, 1000);
	while(gzgets(fi, string, 1000) != NULL){
		sscanf(string,"%s %d %s %s %*s %*s %f %f %*s %s %*s %*s %*s",chr,&Eloc_tmp,ref_kmer,read_id,&Emean_tmp,&Estd_tmp,model_kmer);

		if(strcmp(read_id, last_read_id) == 0 || strcmp(last_read_id, "None")==0){
			Eloc[e_total]=Eloc_tmp;
			Emean[e_total]=Emean_tmp;
			Estd[e_total]=Estd_tmp;
		}else{
			start_site = Eloc[0]; end_site = Eloc[e_total-1]+5;
			cref_seq = faidx_fetch_seq(fai, chr , start_site , end_site , &fetch_len);					
			if(pos_c>neg_c){
				rev_chain=0; rev_mark='+';
			}else{
				rev_chain=1; rev_mark='-';
			}
			cg_region_calling(argv[1], chr, last_read_id, cref_seq, Eloc, Emean, Estd, e_total, start_site, rev_chain, rev_mark);
			
			e_total=0; pos_c=0; neg_c=0;

			Eloc[0]=Eloc_tmp;
			Emean[0]=Emean_tmp;
			Estd[0]=Estd_tmp;
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



	start_site = Eloc[0]; end_site = Eloc[e_total-1]+5;
	cref_seq = faidx_fetch_seq(fai, chr , start_site , end_site , &fetch_len);
	if(pos_c>neg_c){
		rev_chain=0; rev_mark='+';
	}else{
		rev_chain=1; rev_mark='-';
	}
	cg_region_calling(argv[1], chr, read_id, cref_seq, Eloc, Emean, Estd, e_total, start_site, rev_chain, rev_mark);

	return 0;
}


int cg_region_calling(char model_dir[], char chr[], char read_id[], char cref_seq[], int Eloc[], float Emean[], float Estd[], int event_num, int start_site, int rev_chain, char rev_mark){
	int len=strlen(cref_seq);

	char model_p[200],model_a[200],model_cg[200],model_acg[200];
	strcpy(model_p,model_dir); strcat(model_p,"/PCR/PCR.param.txt");
	strcpy(model_a,model_dir); strcat(model_a,"/6mA/6mA.param.txt");
	strcpy(model_cg,model_dir); strcat(model_cg,"/mCG/mCG.param.txt");
	strcpy(model_acg,model_dir); strcat(model_acg,"/6mA_mCG/6mA_mCG.param.txt");
	
	float unmethylscore,cgmethylscore,amethylscore1,amethylscore2,acgmethylscore1,acgmethylscore2;
	int j=0;
	int CG_loc[10000];
	for(long i=0;i<len-1;i++){
		if(cref_seq[i]=='C' && cref_seq[i+1]=='G'){
			CG_loc[j]=i; j++;
		}
	}
	int r_start=0; int r_end=0; int e_start; int e_end; int r_start_t=0; int r_end_t=0;
	char cg_region[1500]; int cg_num=0;
	int shift_pos;
	for(int i=0; i<j; i++){
		if(i==0 && CG_loc[i]>10){
			r_start=CG_loc[i]-10; cg_num=0;
		}
		if(i>0 && CG_loc[i]-CG_loc[i-1]>10){
			if(r_start>0){
				r_end=CG_loc[i-1]+10;
				if(rev_chain==1){
					r_start++;
					r_end++;
					shift_pos=9;
				}else{
					shift_pos=10;
				}
				r_region_to_e_region(Eloc, event_num, r_start,r_end-5, &e_start, &e_end, start_site, &r_start_t, &r_end_t);
				strncpy(cg_region,&cref_seq[r_start_t],r_end_t-r_start_t+1);
				cg_region[r_end_t-r_start_t+1]='\0';

				int e_loc[e_end-e_start+1];
				float e_mean[e_end-e_start+1],e_std[e_end-e_start+1];
				for(int k=0;k<=e_end-e_start;k++)
				{
					e_loc[k]=Eloc[e_start+k];
					e_mean[k]=Emean[e_start+k];
					e_std[k]=Estd[e_start+k];
				}

				unmethylscore=viterbi_score(model_p,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				cgmethylscore=viterbi_score(model_cg,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				amethylscore1=viterbi_score(model_a,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				amethylscore2=viterbi_score(model_a,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,1);
				acgmethylscore1=viterbi_score(model_acg,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				acgmethylscore2=viterbi_score(model_acg,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,1);
				float amethylscore=amethylscore1;
				if(amethylscore<amethylscore2)
					amethylscore=amethylscore2;
				float acgmethylscore=acgmethylscore1;
				if(acgmethylscore<acgmethylscore2)
					acgmethylscore=acgmethylscore2;
				
				printf("%s\t%c\t%d\t%d\t%s\t%.2f\t%s\n",chr,rev_mark,start_site+r_start,start_site+r_end,read_id,acgmethylscore-amethylscore+cgmethylscore-unmethylscore,cg_region);
				//printf("%s\t%c\t%d\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",chr,rev_mark,start_site+r_start,start_site+r_end,read_id,acgmethylscore-cgmethylscore+amethylscore-unmethylscore,acgmethylscore-amethylscore+cgmethylscore-unmethylscore,acgmethylscore1,acgmethylscore2,amethylscore1,amethylscore2,cgmethylscore,unmethylscore,cg_region);
			}
			r_start=CG_loc[i]-10;
			cg_num=0;
		}
		cg_num++;
		if(i==j-1){
			if(len-CG_loc[i]>10){
				r_end=CG_loc[i]+10;
				if(rev_chain==1){
					r_start++;
					r_end++;
					shift_pos=9;
				}else{
					shift_pos=10;
				}
				r_region_to_e_region(Eloc, event_num, r_start,r_end-5, &e_start, &e_end, start_site, &r_start_t, &r_end_t);
				strncpy(cg_region,&cref_seq[r_start_t],r_end_t-r_start_t+1);
				cg_region[r_end_t-r_start_t+1]='\0';
				
				int e_loc[e_end-e_start+1];
				float e_mean[e_end-e_start+1],e_std[e_end-e_start+1];
				for(int k=0;k<=e_end-e_start;k++)
				{
					e_loc[k]=Eloc[e_start+k];
					e_mean[k]=Emean[e_start+k];
					e_std[k]=Estd[e_start+k];
				}
				unmethylscore=viterbi_score(model_p,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				cgmethylscore=viterbi_score(model_cg,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				amethylscore1=viterbi_score(model_a,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				amethylscore2=viterbi_score(model_a,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,1);
				acgmethylscore1=viterbi_score(model_acg,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,0);
				acgmethylscore2=viterbi_score(model_acg,cg_region,e_loc,e_mean,e_std,chr,start_site+r_start_t,e_end-e_start+1,rev_chain,1);
				float amethylscore=amethylscore1;
				if(amethylscore<amethylscore2)
					amethylscore=amethylscore2;
				float acgmethylscore=acgmethylscore1;
				if(acgmethylscore<acgmethylscore2)
					acgmethylscore=acgmethylscore2;


				printf("%s\t%c\t%d\t%d\t%s\t%.2f\t%s\n",chr,rev_mark,start_site+r_start,start_site+r_end,read_id,acgmethylscore-amethylscore+cgmethylscore-unmethylscore,cg_region);
				//printf("%s\t%c\t%d\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",chr,rev_mark,start_site+r_start,start_site+r_end,read_id,acgmethylscore-cgmethylscore+amethylscore-unmethylscore,acgmethylscore-amethylscore+cgmethylscore-unmethylscore,acgmethylscore1,acgmethylscore2,amethylscore1,amethylscore2,cgmethylscore,unmethylscore,cg_region);
			}
		}
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
