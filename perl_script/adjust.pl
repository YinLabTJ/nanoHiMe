use FindBin; 
use File::Basename; 

$perl_main="$FindBin::Bin/$FindBin::Script"; 
my ($filename, $dir) = fileparse($perl_main);

$abd_file=shift; $kmerdir=shift; $outdir=shift;
open IN,"$abd_file";
while(<IN>){
	chomp;
	@t=split;
	$a{$t[0]}=$t[1];
	$b{$t[0]}=$t[2];
}
close IN;


print "$dir/6mer.list\n";
open IN,"$dir/6mer.list";
`mkdir -p $outdir`;
while(<IN>){
	chomp;
	$kmer=$_;
	open OUT,">$outdir/$kmer.eventalign.txt";
	open EVENT,"$kmerdir/$kmer.eventalign.txt";
	while(<EVENT>){
		chomp;
		@t=split;
		if(defined($b{$t[0]}) && $b{$t[0]} ne "NA"){
			$mean=($t[1]-$a{$t[0]})/$b{$t[0]};
			print OUT "$t[0]\t$mean\n" if($t[3]>0.00149);
		}
	}
	close EVENT;
	close OUT;
}
close IN;
