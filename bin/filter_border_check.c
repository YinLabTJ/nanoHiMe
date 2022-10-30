#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{
	FILE *ep;
	char command1[100];
	//sprintf(command1,"mkdir -p 6mers_c/%s",argv[2]);
	//system(command1);
	if((ep = fopen(argv[1],"r")) == NULL)
	{
		printf("Can't open event file %s\n",argv[3]);
		exit(EXIT_FAILURE);
	}
	
	int Eloc_tmp; int last_loc=0;
	int e_total=0;
	float Emean[100000],Estd[100000]; int Eloc[100000]; float Eduration[100000];
	char K[100000][7];
	float Emean_tmp,Estd_tmp,Eduration_tmp;
	char read_id[50]; char last_read_id[50]="None";
	char kmer[7];
	int Readnum[100000]; int readn;

	while(fscanf(ep,"%*s %d %*s %s %*s %d %f %f %f %s %*s %*s %*s",&Eloc_tmp,read_id,&readn,&Emean_tmp,&Estd_tmp,&Eduration_tmp,kmer)!=EOF){
		//printf("%d\t%s\t%f\t%f\t%f\t%s\n",Eloc_tmp,read_id,Emean_tmp,Estd_tmp,Eduration_tmp,kmer);
		if((strcmp(read_id, last_read_id) == 0)|| strcmp(last_read_id, "None")==0){
			Eloc[e_total]=Eloc_tmp;
			last_loc=Eloc_tmp;
			Readnum[e_total]=readn;
			Emean[e_total]=Emean_tmp;
			Estd[e_total]=Estd_tmp;
			Eduration[e_total]=Eduration_tmp;
			strcpy(K[e_total], kmer);
		}else{
			if(e_total>10){
				for(int i=5;i<=e_total-5;i++){
					//char command_inf[200];
					printf("%s\t%s\t%.2f\t%.2f\t%f\t%d\n",last_read_id,K[i],Emean[i],Estd[i],Eduration[i],Readnum[i]);
					//sprintf(command_inf,"echo %s %.2f %.2f %f >> 6mers_c/%s/%s.eventalign.txt",last_read_id,Emean[i],Estd[i],Eduration[i],argv[2],K[i]);
					//system(command_inf);
				}
			}
			e_total=0; 
			Eloc[0]=Eloc_tmp;
			last_loc=Eloc_tmp;
			Emean[0]=Emean_tmp;
			Estd[0]=Estd_tmp;
			Eduration[0]=Eduration_tmp;
			strcpy(K[0], kmer);
		}
		e_total++;
		strcpy(last_read_id, read_id);

	}
	fclose(ep);

	if(e_total>10){
		for(int i=5;i<=e_total-5;i++){
			printf("%s\t%s\t%.2f\t%.2f\t%f\t%d\n",last_read_id,K[i],Emean[i],Estd[i],Eduration[i],Readnum[i]);
			//char command_inf[200];
			//sprintf(command_inf,"echo %s %.2f %.2f %f >> 6mers_c/%s/%s.eventalign.txt",last_read_id,Emean[i],Estd[i],Eduration[i],argv[2],K[i]);
			//system(command_inf);
		}
	}

	return 0;
}
