open IN,"PCR/PCR.param.txt";
while(<IN>){
	chomp;
	@t=split;
	$pcr{$t[0]}=$t[1];
}
close IN;
open IN,"mCG/mCG.param.txt";
while(<IN>){
	chomp;
	next if(/\#/);
	@t=split;
	$mcg{$t[0]}=$t[1];
}
close IN;
open OUT,">PCR_vs_mCG.input";
print OUT "Kmer\tPCR\tmCG\tClass\n";
foreach $key(sort keys %pcr){
	if($key=~/CG/){
		$class="CG";
	}elsif($key=~/C$/){
		$class="XXXXXC";
	}else{
		$class="nonCG";
	}
	print OUT "$key\t$pcr{$key}\t$mcg{$key}\t$class\n";
}
close OUT;
