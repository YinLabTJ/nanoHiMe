use FindBin;
use File::Basename;

my $perl_main="$FindBin::Bin/$FindBin::Script";
my ($filename, $dir) = fileparse($perl_main);

$mean_dir=shift; $final_kmer_dir=shift; $old_model=shift;
open IN,"$old_model";
while(<IN>){
	chomp;
	@t=split;
	$old_mean{$t[0]}=$t[1];
	$old_stdv{$t[0]}=$t[2];
}
close IN;

open IN,"$mean_dir/filtered.num.txt";
while(<IN>){
	chomp;
	@t=split /\//;
	$t[-1]=~s/\.txt//g;
	$num=<IN>; chomp($num);
	$filter{$t[-1]}=1 if($num<30);
}
close IN;


open KMER,"$dir/6mer.list";
while(<KMER>){
	chomp;
	$kmer=$_;
	open IN,"$final_kmer_dir/$kmer.eventalign.txt";
	$i=0; $sum=0; $n2=0; @E="";
	while(<IN>){
		chomp;
		@t=split;
		next if($filter{$t[0]});
		$E[$i]=$t[1];
		$sum+=$E[$i];
		$i++;
	}
	close IN;
	if($i>50){
		$mean=$sum/$i; 
		for($j=0;$j<$i;$j++){
			$n2+=(($E[$j]-$mean)*($E[$j]-$mean));
		}
		$stdv=$n2/($i-1);
		$stdv=$stdv**0.5;
		$mean=sprintf "%.2f",$mean;
		$stdv=sprintf "%.2f",$stdv;
	}else{
		$mean=$old_mean{$kmer};
		$stdv=$old_stdv{$kmer}
	}
	print "$kmer\t$mean\t$stdv\n";
}
close KMER;
