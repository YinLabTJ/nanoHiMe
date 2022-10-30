open IN,"6mA.double.xls";
while(<IN>){
	chomp;
	@t=split;
	$filter{$t[0]}=1;
}
close IN;
open IN,"6mA_mCG.double.xls";
while(<IN>){
	chomp;
	@t=split;
	$filter{$t[0]}=1;
}
close IN;

open IN,"6mA/6mA.param.txt";
while(<IN>){
	chomp;
	@t=split;
	next if($filter{$t[0]});
	$ma6{$t[0]}=$t[1];
}
close IN;
open IN,"6mA_mCG/6mA_mCG.param.txt";
while(<IN>){
	chomp;
	next if(/\#/);
	@t=split;
	next if($filter{$t[0]});
	$ma6_mcg{$t[0]}=$t[1];
}
close IN;
open OUT,">6mA_vs_6mA_mCG.input";
print OUT "Kmer\t6mA\t6mA_mCG\tClass\n";
foreach $key(sort keys %ma6){
	if($key=~/CG/){
		$class="CG";
	}elsif($key=~/C$/){
		$class="XXXXXC";
	}else{
		$class="nonCG";
	}
	print OUT "$key\t$ma6{$key}\t$ma6_mcg{$key}\t$class\n";
}
close OUT;
