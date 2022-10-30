#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

int main(int argc, char *argv[])
{
	/*FILE *ep;
	if((ep = fopen(argv[1],"r")) == NULL)
	{
		printf("Can't open event file %s\n",argv[1]);
		exit(EXIT_FAILURE);
	} */
	
	int Eloc_tmp; int last_loc=0;
	int number=0;
	float Emean_tmp,Estd_tmp,Eduration_tmp;
	char read_id[50],chr[50]; char last_read_id[50]="None";
	char kmer[7];

	gzFile fi;
	fi = gzopen(argv[1],"rb");
	char string[1000];

	gzgets(fi, string, 1000);
	while(gzgets(fi, string, 1000) != NULL){
		sscanf(string,"%s %d %*s %s %*s %*s %f %f %f %s %*s %*s %*s",chr,&Eloc_tmp,read_id,&Emean_tmp,&Estd_tmp,&Eduration_tmp,kmer);
		//printf("%s",string);
		if(strcmp(last_read_id, "None")==0){
			printf("%s\t%d\t",chr,Eloc_tmp);
		}else{
			if(strcmp(read_id, last_read_id) == 0){
			}else{
				number++;
				printf("%d\t%s\t%d\n%s\t%d\t",last_loc,last_read_id,number,chr,Eloc_tmp);
			}
		}
		last_loc=Eloc_tmp;
		strcpy(last_read_id, read_id);

	}
	
	number++;
	printf("%d\t%s\t%d\n",last_loc,last_read_id,number);

	/*fscanf(ep,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	while(fscanf(ep,"%s %d %*s %s %*s %*s %f %f %f %s %*s %*s %*s",chr,&Eloc_tmp,read_id,&Emean_tmp,&Estd_tmp,&Eduration_tmp,kmer)!=EOF){
		//printf("%d\t%s\t%f\t%f\t%f\t%s\n",Eloc_tmp,read_id,Emean_tmp,Estd_tmp,Eduration_tmp,kmer);
		//if((strcmp(read_id, last_read_id) == 0)|| strcmp(last_read_id, "None")==0){
			if(strcmp(last_read_id, "None")==0){
				printf("%s\t%d\t",chr,Eloc_tmp);
			}else{
				if(strcmp(read_id, last_read_id) == 0){
				}else{
					number++;
					printf("%d\t%s\t%d\n%s\t%d\t",last_loc,last_read_id,number,chr,Eloc_tmp);
				}
			}
			last_loc=Eloc_tmp;
			strcpy(last_read_id, read_id);
	}
	fclose(ep);

	number++;
	printf("%d\t%s\t%d\n",last_loc,last_read_id,number); */
}
