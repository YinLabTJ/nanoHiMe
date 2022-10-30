#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "kmer2binary.h"

#define pi 3.141592653589

struct Kmer_model
{
	char kmer[7];
	float mean;
	float std;
};

float viterbi_score(char model_file[], char input_sequence[], int Eloc[], float Emean[], float Estd[], char chr[], int start_site, int event_total, int rev_chain, int methyl_state, char read_id[], float Eduration[])
{
	int ch;
	struct Kmer_model PCR[4096];
	FILE *mp; //model_file
	FILE *ep; //event_file
	FILE *op; //output_file

	float final_score=0;

		char ref_sequence[30000];
		if(methyl_state==1){
			if(rev_chain==0)
				methylation_transition(input_sequence, ref_sequence);
			if(rev_chain==1)
				methylation_transition_neg(input_sequence, ref_sequence);
		}else{
			strcpy(ref_sequence,input_sequence);
		}
		//puts(ref_sequence);
		if((mp = fopen(model_file,"r")) == NULL)
		{
			printf("Can't open %s\n", model_file);
			exit(EXIT_FAILURE);
		}
		for(int i=0; i < 4096; i++){
			fscanf(mp,"%s %f %f",PCR[i].kmer, &PCR[i].mean, &PCR[i].std);
			//printf("%d %s %f %f\n",i, PCR[i].kmer, PCR[i].mean, PCR[i].std);
		}
		fclose(mp);
		//printf("Model %s has been loaded\n", model_file);
			


		//ref_seq
		int len = strlen(ref_sequence);
		//long event_total = atoi(argv[6]);
		float Rmean[len-5],Rstd[len-5]; //float Emean[event_total],Estd[event_total];
		//long Eloc[event_total];

		//printf("The length of Sequence %s is %d.\n", ref_sequence, len);
		char seq1[7]; char seq2[7]; char seq1b[8]; char seq2b[8];
		for(int i=1; i<=len-5; i++){
			strncpy(seq1,&ref_sequence[i-1],6);
			seq1[6]='\0';
			float CG_bias;

			if(rev_chain==1){
				reverse_complement(seq1, seq2);
				/*if(methyl_state==1){
					char seqt[7];
					strcpy(seqt,seq2);
					methylation_transition(seqt,seq2);
				}*/
				int num = kmer_to_num(seq2);

				if(i>1){
					strncpy(seq1b,&ref_sequence[i-2],7);
					seq1b[7]='\0';
					reverse_complement(seq1b, seq2b);
					CG_bias=includeCG(seq2b);
				}else{
					CG_bias=0;
				}
				Rmean[i-1]=PCR[num].mean;
				Rstd[i-1]=PCR[num].std+CG_bias;
			}else{
				int num = kmer_to_num(seq1);
				if(i<len-5){
					strncpy(seq1b,&ref_sequence[i-1],7);
					CG_bias=includeCG(seq1b);
				}else{
					CG_bias=0;
				}
				Rmean[i-1]=PCR[num].mean;
				Rstd[i-1]=PCR[num].std+CG_bias;
			}
		}

		//printf("Matrice has %d 6-mers and %d events\n",len-5,event_total);	
	


		/*int Event_round_max=50;
		int Ref_round_max=100;

		char ***S_stat;
		int si,sj;
		S_stat=(char ***)malloc(Event_round_max * sizeof(char **));
		for(sj=0;sj<Event_round_max;++sj)
			S_stat[sj]=(char **)malloc(Event_round_max*3 * sizeof(char *));
		//printf("%d * %d * %d\n",Event_round_max,Ref_round_max,event_total+len);
		for(sj=0;sj<Event_round_max;++sj){
			for(si=0;si<Ref_round_max*3;si++){
				S_stat[sj][si]=(char *)malloc((event_total+len)*2*sizeof(char));
			}
		}
		*/	

		

		//printf("Memory prepared\n");
		

		char final_path[(event_total + len)*2];
		final_path[0]='\0';

		//float S_loglik[Event_round_max][Ref_round_max*3];
		//float S_loglik[event_total][(len-5)*3];
		float S_loglik[50][300];

		//char S_stat[event_total][(len-5)*3][(event_total+len)*2];
		char S_stat[50][300][150];

		int e_start=0; int r_start=0; int max_col; int start_state=2;
		

		while(e_start<event_total-1){
				int e_end = e_start+49 >= event_total-1 ? event_total-1 : e_start+49;
				int r_end;
				
				r_end = Eloc[e_end]-start_site;
				if(r_end <= Eloc[event_total-1]-start_site-2)
					r_end=r_end+2;
				//r_end = r_start + 99 > Eloc[event_total-1]-start_site ? Eloc[event_total-1]-start_site : r_start+99;

				//printf("Event %d to %d and Ref %d to %d      %d - %d\n",e_start,e_end,r_start,r_end,Eloc[e_end],start_site);
				for(int i=0;i<=r_end-r_start;i++){
					for(int j=0;j<=e_end-e_start;j++){
						S_loglik[j][i*3]=-INFINITY; S_stat[j][i*3][0]='\0';
						S_loglik[j][i*3+1]=-INFINITY; S_stat[j][i*3+1][0]='\0';
						S_loglik[j][i*3+2]=-INFINITY; S_stat[j][i*3+2][0]='\0';
					}
				}

				//printf("Initialing finished.\n");

				//printf("start_stat = %d\n",start_state);
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
				//printf("Viterbi Matrix finish\n");
		
				/*for(int i=0;i<=1;i++){
					for(int j=0;j<=1;j++){
						printf("no %d event, no %d 6mer is \n %f and %s \n %f and  %s \n %f and %s \n",i+1,j+1,S_loglik[i][j*3], S_stat[i][j*3], S_loglik[i][j*3+1], S_stat[i][j*3+1], S_loglik[i][j*3+2], S_stat[i][j*3+2]);
					}	
				}*/
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
					//if(atoi(&S_stat[e_end-e_start][max_col][path_len-1-end_skip])>=7){
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
						//start_state = atoi(&S_stat[e_end-e_start][max_col][path_len-1-end_skip])<=4 ? 2 : 1;
					}
				}
				//printf("%s , skip %d and state %d\n",S_stat[e_end-e_start][max_col],end_skip,start_state);
				
				S_stat[e_end-e_start][max_col][path_len-end_skip]='\0';
				
				//printf("final_path0: %s\n",final_path);
				strcat(final_path,S_stat[e_end-e_start][max_col]);
				//printf("final_path: %s\n",final_path);
				if(end_skip==0){
					final_score+=S_loglik[e_end-e_start][max_col];
				}else{
					final_score+=S_loglik[e_end-e_start][col1-end_skip*3+start_state];
				}
				long fpath_len=strlen(final_path);

				//printf("e_start: %d, e_end: %d , final_score : %.2f, end_skip: %d \n",e_start,e_end,final_score,end_skip);

				e_start+=49;
		
		}
		
		
		char rev_sequence[len+1];
		reverse_complement(ref_sequence, rev_sequence);
		
		int j=0; int k=0;
		
		//char output_info[event_total][150];

		if(rev_chain==1){
			//strncpy(seq1,&ref_sequence[j],6);
			//seq1[6]='\0';
			//reverse_complement(seq1,seq2);
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

		//if(strcmp(input_sequence,"ATATGCACACGTATGTTTATT")==0)
			//printf("%s\t%d\t%s\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\n",chr,start_site,seq1,Emean[0],Estd[0],seq2,PCR[num].mean,PCR[num].std);

		printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq1,read_id,0,Emean[0],Estd[0],Eduration[0],seq2,PCR[num].mean,PCR[num].std);

		//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\t%s\n",chr,start_site,Emean[0],Estd[0],seq2,PCR[num].mean,PCR[num].std,seq1);
		long path_len=strlen(final_path);
		//printf("final path is %s\nlength is %ld\n",final_path,path_len);
		for(long i=0;i<path_len;i++){
			if(final_path[i]=='0'){
				k=k+1;
				if(rev_chain==1){
					//strncpy(seq1,&rev_sequence[j],6);
					//seq1[6]='\0';
					//reverse_complement(seq1,seq2);
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
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\t%s\n",chr,start_site,Emean[k],Estd[k],seq2,PCR[num].mean,PCR[num].std,seq1);
			}
			if(final_path[i]=='1'){
				j=j+1; k=k+1;
				start_site=start_site+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6);
					seq1[6]='\0';
					strncpy(seq2,&rev_sequence[len-j-6],6);
					seq2[6]='\0';
					//strncpy(seq1,&rev_sequence[j],6);
					//seq1[6]='\0';
					//reverse_complement(seq1,seq2);
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq2);
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);	
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\t%s\n",chr,start_site,Emean[k],Estd[k],seq2,PCR[num].mean,PCR[num].std,seq1);
			}
			if(final_path[i]=='2'){
				k=k+1;
				if(rev_chain==1){
					//strncpy(seq1,&rev_sequence[j],6);
					//seq1[6]='\0';
					//reverse_complement(seq1,seq2);
					strncpy(seq1,&input_sequence[j],6); seq1[6]='\0'; strncpy(seq2,&rev_sequence[len-j-6],6); seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq2);
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\t%s\n",chr,start_site,Emean[k],Estd[k],seq2,PCR[num].mean,PCR[num].std,seq1);
			}
			if(final_path[i]=='3'){
				j=j+1; k=k+1;
				start_site=start_site+1;
				if(rev_chain==1){
					/*strncpy(seq1,&rev_sequence[j],6);
					seq1[6]='\0';
					reverse_complement(seq1,seq2);*/
					strncpy(seq1,&input_sequence[j],6); seq1[6]='\0'; strncpy(seq2,&rev_sequence[len-j-6],6); seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}	
				int num = kmer_to_num(seq2);
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\t%s\n",chr,start_site,Emean[k],Estd[k],seq2,PCR[num].mean,PCR[num].std,seq1);
			}
			if(final_path[i]=='4'){
				j=j+1; k=k+1;
				start_site=start_site+1;
				if(rev_chain==1){
					/*strncpy(seq1,&rev_sequence[j],6);
					seq1[6]='\0';
					reverse_complement(seq1,seq2);*/
					strncpy(seq1,&input_sequence[j],6); seq1[6]='\0'; strncpy(seq2,&rev_sequence[len-j-6],6); seq2[6]='\0';
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);
				}
				int num = kmer_to_num(seq2);
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);	
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t%3.2f\t%3.2f\t%s\n",chr,start_site,Emean[k],Estd[k],seq2,PCR[num].mean,PCR[num].std,seq1);
			}
			if(final_path[i]=='5'){
				k=k+1;
				if(rev_chain==1){
					//strncpy(seq1,&rev_sequence[j],6);
					strncpy(seq1,&input_sequence[j],6); 
				}else{
					strncpy(seq1,&ref_sequence[j],6);
				}
				seq1[6]='\0';
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\tNNNNNN\t0.00\t0.00\tinf\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k]);
					//printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq2,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t0\t0\tNNNNNN\n",chr,start_site,Emean[k],Estd[k],seq2);
			}
			if(final_path[i]=='6'){
				k=k+1;
				if(rev_chain==1){
					strncpy(seq1,&input_sequence[j],6);
					/*strncpy(seq1,&rev_sequence[j],6);
					seq1[6]='\0';
					reverse_complement(seq1,seq2);*/
				}else{
					strncpy(seq1,&ref_sequence[j],6);
					/*strncpy(seq1,&ref_sequence[j],6);
					seq1[6]='\0';
					strcpy(seq2,seq1);*/
				}
				seq1[6]='\0';
					printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\tNNNNNN\t0.00\t0.00\tinf\n",chr,start_site,seq1,read_id,k,Emean[k],Estd[k],Eduration[k]);
					//printf("%s\t%d\t%s\t%s\tt\t%d\t%3.2f\t%1.3f\t%1.5f\t%s\t%3.2f\t%3.2f\t0\n",chr,start_site,seq2,read_id,k,Emean[k],Estd[k],Eduration[k],seq2,PCR[num].mean,PCR[num].std);
				//sprintf(output_info[k],"%s\t%d\t%3.2f\t%3.2f\t%s\t0\t0\tNNNNNN\n",chr,start_site,Emean[k],Estd[k],seq2);
			}
			if(final_path[i]=='7'){
				//if(rev_chain==1){
				//	start_site=start_site-1;
				//}else{
					start_site=start_site+1;
				//}
				j=j+1;
			}
			if(final_path[i]=='8'){
				//if(rev_chain==1){
				//	start_site=start_site-1;
				//}else{
					start_site=start_site+1; 
				//}
				j=j+1;
			}
			if(final_path[i]=='9'){
				//if(rev_chain==1){
				//	start_site=start_site-1;
				//}else{
					start_site=start_site+1;
				//}
				j=j+1;
			}
		}
		
		//printf("Free Memory %d %d\n",Event_round_max,Ref_round_max*3);
		//for(sj=0;sj<Event_round_max;++sj){
		//	for(si=0;si<Ref_round_max*3;si++){
		//		printf("sj:%d si:%d\n",sj,si);
		//		free((void *)S_stat[sj][si]);
                //        }
		//}
		//for(sj=0;sj<Event_round_max;++sj)
		//	free(S_stat[sj]);
		//free(S_stat);
		
		
		//if((op=fopen("our_output.txt","w"))==NULL){
		//	fprintf(stderr,"Can't create ouput file.\n");
		//}
		//for(long i=0; i<event_total; i++)
			//fprintf(op,"%s",output_info[i]);
	return final_score;
}


