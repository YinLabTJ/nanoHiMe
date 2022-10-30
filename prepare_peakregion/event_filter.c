#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>


int main(int argc, char *argv[])
{
	gzFile fi;
	fi = gzopen(argv[1],"rb");
	
	FILE *op;
	if((op = fopen(argv[2],"r")) == NULL)
	{
		printf("Can't open overlap file %s\n",argv[2]);
		exit(EXIT_FAILURE);
	}


	int start_tmp,end_tmp,read_number;
	int start[1000000]={0};
	int end[1000000]={0};
	while(fscanf(op,"%*s %d %d %*s %d",&start_tmp,&end_tmp,&read_number)!=EOF){
		if(start[read_number]==0 || start[read_number]>start_tmp)
			start[read_number]=start_tmp;
		if(end[read_number]==0 || end[read_number]<end_tmp)
			end[read_number]=end_tmp;
	}
	fclose(op);

	int Eloc_tmp,Enum_tmp; int last_loc=0;
	int number=0;
	char Emean_tmp[15],Estd_tmp[15],Eduration_tmp[15];
	char read_id[50],chr[50]; char last_read_id[50]="None";
	char model_kmer[7],kmer[7];
	char v1[10],v2[10],v3[10];

	char string[1000];

	gzgets(fi, string, 1000);
	while(gzgets(fi, string, 1000) != NULL){
		sscanf(string,"%s %d %s %s %*s %d %s %s %s %s %s %s %s",chr,&Eloc_tmp,model_kmer,read_id,&Enum_tmp,Emean_tmp,Estd_tmp,Eduration_tmp,kmer,v1,v2,v3);
		if(strcmp(last_read_id, "None")==0){
			number=1;
		}else{
			if(strcmp(read_id, last_read_id) != 0){
				number++;
			}
		}
		if(start[number]>0){
			if(Eloc_tmp>=start[number] && Eloc_tmp<=end[number]){
				printf("%s\t%d\t%s\t%s\tt\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",chr,Eloc_tmp,model_kmer,read_id,Enum_tmp,Emean_tmp,Estd_tmp,Eduration_tmp,kmer,v1,v2,v3);
			}
		}
		strcpy(last_read_id, read_id);
	}

}
