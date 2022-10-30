#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "kmer2binary_cg.h"

#define pi 3.141592653589

struct Kmer_model
{
	char kmer[7];
	float mean;
	float std;
	float mean2;
	float std2;
};

float viterbi_score(char model_file[], char input_sequence[], int Eloc[], float Emean[], float Estd[], char chr[], int start_site, int event_total, int rev_chain, int methyl_state)
{
	int ch;
	struct Kmer_model PCR[4096];
	FILE *mp; //model_file
	FILE *ep; //event_file
	FILE *op; //output_file

	float final_score=0;

		char ref_sequence[30000];
		strcpy(ref_sequence, input_sequence);
		
		if((mp = fopen(model_file,"r")) == NULL)
		{
			printf("Can't open %s\n", model_file);
			exit(EXIT_FAILURE);
		}
		for(int i=0; i < 4096; i++){
			fscanf(mp,"%s %f %f %f %f",PCR[i].kmer, &PCR[i].mean, &PCR[i].std, &PCR[i].mean2, &PCR[i].std2);
		}
		fclose(mp);
			


		//ref_seq
		int len = strlen(ref_sequence);
		float Rmean[len-5],Rstd[len-5]; 

		char seq1[7]; char seq2[7];
		for(int i=1; i<=len-5; i++){
			strncpy(seq1,&ref_sequence[i-1],6);
			seq1[6]='\0';
			if(rev_chain==1){
				reverse_complement(seq1, seq2);
				int num = kmer_to_num(seq2);
				if(methyl_state==0){
					Rmean[i-1]=PCR[num].mean;
					Rstd[i-1]=PCR[num].std;
				}else{
					Rmean[i-1]=PCR[num].mean2;
					Rstd[i-1]=PCR[num].std2;
				}
			}else{
				int num = kmer_to_num(seq1);
				if(methyl_state==0){
					Rmean[i-1]=PCR[num].mean;
					Rstd[i-1]=PCR[num].std;
				}else{
					Rmean[i-1]=PCR[num].mean2;
					Rstd[i-1]=PCR[num].std2;
				}
			}
		}


	


		char final_path[(event_total + len)*2];
		final_path[0]='\0';

		float S_loglik[50][150];

		char S_stat[50][150][120];

		int e_start=0; int r_start=0; int max_col; int start_state=2;
		

		while(e_start<event_total-1){
				int e_end = e_start+49 >= event_total-1 ? event_total-1 : e_start+49;
				int r_end;
				
				r_end = Eloc[e_end]-start_site;
				
				for(int i=0;i<=r_end-r_start;i++){
					for(int j=0;j<=e_end-e_start;j++){
						S_loglik[j][i*3]=-INFINITY; S_stat[j][i*3][0]='\0';
						S_loglik[j][i*3+1]=-INFINITY; S_stat[j][i*3+1][0]='\0';
						S_loglik[j][i*3+2]=-INFINITY; S_stat[j][i*3+2][0]='\0';
					}
				}


				if(start_state==2){
					S_loglik[0][2]=-1*(pow(Emean[e_start]-Rmean[r_start],2))/(2*pow(Rstd[r_start],2))+log(1/pow((2*pi*pow(Rstd[r_start],2)),0.5));
				}else{
					S_loglik[0][1]=0;
				}

				for(int i=0;i<=r_end-r_start;i++){
					for(int j=1;j<=e_end-e_start;j++){
						
						//Match_self
						char movement[2];
						int t1=r_end-r_start+1; int t2=e_end-e_start+1;
						float events_per_base=t2/t1;
						if(events_per_base<1.25)
							events_per_base=1.25;
						float S0=S_loglik[j-1][i*3+2]+log(1-1/events_per_base)+(-1*(pow((Emean[e_start+j]-Rmean[r_start+i]),2))/(2*pow(Rstd[r_start+i],2)))+log(1/pow(2*pi*pow(Rstd[r_start+i],2),0.5));
						S_loglik[j][i*3+2]=S0; strcpy(S_stat[j][i*3+2],S_stat[j-1][i*3+2]);  strcpy(movement,"0");
						
						//Match prev
						if(i>0){
							float S1=S_loglik[j-1][i*3-1]+log(1-(1-1/events_per_base)-0.001-0.0025)+(-1*(pow((Emean[e_start+j]-Rmean[r_start+i]),2))/(2*pow(Rstd[r_start+i],2)))+log(1/pow(2*pi*pow(Rstd[r_start+i],2),0.5));
							float tmp1=log(1-(1-(len-5)/event_total)-0.001-0.0025);
							if(S1>S_loglik[j][i*3+2]){
								S_loglik[j][i*3+2]=S1; strcpy(S_stat[j][i*3+2],S_stat[j-1][i*3-1]); strcpy(movement,"1");
							}
						}
						
						//Bad_self
						float S2=S_loglik[j-1][i*3+1]+log(0.333)+(-1*(pow((Emean[e_start+j]-Rmean[r_start+i]),2))/(2*pow(Rstd[r_start+i],2)))+log(1/pow(2*pi*pow(Rstd[r_start+i],2),0.5));
						if(S2>S_loglik[j][i*3+2]){
							S_loglik[j][i*3+2]=S2; strcpy(S_stat[j][i*3+2],S_stat[j-1][i*3+1]); strcpy(movement,"2");
						}
						//Bad_prev
						if(i>0){
							float S3=S_loglik[j-1][i*3-2]+log(0.333)+(-1*(pow((Emean[e_start+j]-Rmean[r_start+i]),2))/(2*pow(Rstd[r_start+i],2)))+log(1/pow(2*pi*pow(Rstd[r_start+i],2),0.5));
							if(S3>S_loglik[j][i*3+2]){
								S_loglik[j][i*3+2]=S3; strcpy(S_stat[j][i*3+2],S_stat[j-1][i*3-2]); strcpy(movement,"3");
							}
						}
						//skip to Match
						if(i>0){
							float S4=S_loglik[j-1][i*3-3]+log(0.7)+(-1*(pow((Emean[e_start+j]-Rmean[r_start+i]),2))/(2*pow(Rstd[r_start+i],2)))+log(1/pow(2*pi*pow(Rstd[r_start+i],2),0.5));
							if(S4>S_loglik[j][i*3+2]){
								S_loglik[j][i*3+2]=S4; strcpy(S_stat[j][i*3+2],S_stat[j-1][i*3-3]); strcpy(movement,"4");
							}
						}
						strcat(S_stat[j][i*3+2],movement);
						//Match_self_Bad
						float S5=S_loglik[j-1][i*3+2]+log(0.001)+0;
						S_loglik[j][i*3+1]=S5; strcpy(S_stat[j][i*3+1],S_stat[j-1][i*3+2]); strcpy(movement,"5");

						//Bad_self_Bad
						float S6=S_loglik[j-1][i*3+1]+log(0.001)+0;
						if(S6>S_loglik[j][i*3+1]){
							S_loglik[j][i*3+1]=S6; strcpy(S_stat[j][i*3+1],S_stat[j-1][i*3+1]); strcpy(movement,"6");
						}
						strcat(S_stat[j][i*3+1],movement);
				

						//Match_prev_Skip
						if(i>0){
							float S7=S_loglik[j][i*3-1]+log(0.0025)+0;
							S_loglik[j][i*3]=S7; strcpy(S_stat[j][i*3],S_stat[j][i*3-1]); strcpy(movement,"7");
						//Bad_prev_Skip
							float S8=S_loglik[j][i*3-2]+log(0.333)+0;

							if(S8>S_loglik[j][i*3]){
								S_loglik[j][i*3]=S8; strcpy(S_stat[j][i*3],S_stat[j][i*3-2]); strcpy(movement,"8");
							}
				

						//Skip prev Skip
				
							float S9=S_loglik[j][i*3-3]+log(0.3)+0;
							if(S9>S_loglik[j][i*3]){
								S_loglik[j][i*3]=S9; strcpy(S_stat[j][i*3],S_stat[j][i*3-3]); strcpy(movement,"9");
							}
							strcat(S_stat[j][i*3],movement);
						}

					}
				}
				
				int col1=(r_end-r_start)*3;
				int col2=(r_end-r_start)*3+1;
				int col3=(r_end-r_start)*3+2;
				
				float max_end = S_loglik[e_end-e_start][col1];
				max_col = col1;
				if(S_loglik[e_end-e_start][col2] > max_end){
					max_end = S_loglik[e_end-e_start][col2];
					max_col = col2;
				}
				if(S_loglik[e_end-e_start][col3] > max_end){
					max_end = S_loglik[e_end-e_start][col3];
					max_col = col3;
				}

				int r_tmp=r_start;	
				r_start=r_end;
				int path_len=strlen(S_stat[e_end-e_start][max_col]);
				int mark=0;
				int end_skip=0;
				while(mark==0){
					if(S_stat[e_end-e_start][max_col][path_len-1-end_skip]=='7' || S_stat[e_end-e_start][max_col][path_len-1-end_skip]=='8' || S_stat[e_end-e_start][max_col][path_len-1-end_skip]=='9'){
						r_start--;
						end_skip++;
					}else{
						mark=1;
						
						if(S_stat[e_end-e_start][max_col][path_len-1-end_skip]=='5' || S_stat[e_end-e_start][max_col][path_len-1-end_skip]=='6'){
							start_state = 1;
						}else{
							start_state = 2;
						}
					}
				}
				
				S_stat[e_end-e_start][max_col][path_len-end_skip]='\0';
				
				strcat(final_path,S_stat[e_end-e_start][max_col]);
				if(end_skip==0){
					final_score+=S_loglik[e_end-e_start][max_col];
				}else{
					final_score+=S_loglik[e_end-e_start][col1-end_skip*3+start_state];
				}
				long fpath_len=strlen(final_path);


				e_start+=49;
		
		}
		
		
		char rev_sequence[len+1];
		reverse_complement(ref_sequence, rev_sequence);
		
		int j=0; int k=0;

		
		char output_info[event_total][150];

		if(rev_chain==1){
			strncpy(seq1,&input_sequence[j],6);
			seq1[6]='\0';
			strncpy(seq2,&rev_sequence[len-j-6],6);
			seq2[6]='\0';

		}else{
			strncpy(seq1,&ref_sequence[j],6);
			seq1[6]='\0';
			strcpy(seq2,seq1);
		}
		int num = kmer_to_num(seq2);

		long path_len=strlen(final_path);
		for(long i=0;i<path_len;i++){
			if(final_path[i]=='0'){
				k=k+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6);
					seq1[6]='\0';
					strncpy(seq2,&rev_sequence[len-j-6],6);
					seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq2);
			}
			if(final_path[i]=='1'){
				j=j+1; k=k+1;
				start_site=start_site+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6);
					seq1[6]='\0';
					strncpy(seq2,&rev_sequence[len-j-6],6);
					seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq2);
			}
			if(final_path[i]=='2'){
				k=k+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6); seq1[6]='\0'; strncpy(seq2,&rev_sequence[len-j-6],6); seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq1);
			}
			if(final_path[i]=='3'){
				j=j+1; k=k+1;
				start_site=start_site+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6); seq1[6]='\0'; strncpy(seq2,&rev_sequence[len-j-6],6); seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}	
				int num = kmer_to_num(seq2);
			}
			if(final_path[i]=='4'){
				j=j+1; k=k+1;
				start_site=start_site+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6); seq1[6]='\0'; strncpy(seq2,&rev_sequence[len-j-6],6); seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq2);
			}
			if(final_path[i]=='5'){
				k=k+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6); 
				}else{
					strncpy(seq1,&ref_sequence[j],6);
				}
				seq1[6]='\0';
			}
			if(final_path[i]=='6'){
				k=k+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6);
				}else{
					strncpy(seq1,&ref_sequence[j],6);
				}
				seq1[6]='\0';
			}
			if(final_path[i]=='7'){
					start_site=start_site+1;
				j=j+1;
			}
			if(final_path[i]=='8'){
				start_site=start_site+1; 
				j=j+1;
			}
			if(final_path[i]=='9'){
				start_site=start_site+1;
				j=j+1;
			}
		}
		
	return final_score;
}


