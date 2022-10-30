struct Kmer_model
{
	char kmer[7];
	float mean;
	float std;
};

//float viterbi_score(char model_file[], char ref_sequence[], char event_file[], char chr[], long start_site, long event_total, int rev_chain, int methyl_state);
float viterbi_score(char model_file[], char ref_sequence[], int Eloc[], float Emean[], float Estd[], char chr[], int start_site, int event_total, int rev_chain, int methyl_state, char read_id[], float Eduration[]);
