$fasta=shift;
open IN,"$fasta";
while(<IN>){
	chomp;
	if(/^>/){
		print "$_\n";
	}else{
		$_=~tr/atcg/ATCG/;
		print "$_\n";
	}
}
close IN; close OUT;
