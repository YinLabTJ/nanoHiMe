$infile=shift; $outdir=shift;
`mkdir -p $outdir`;
open IN,"$infile";
$line=<IN>;
@t=split /\t/,$line;
open OUT,">$outdir/$t[1].eventalign.txt";
$kmer=$t[1];
print OUT "$t[0]\t$t[2]\t$t[3]\t$t[4]\n";
while(<IN>){
	chomp;
	@t=split;
	if($t[1] ne $kmer){
		close OUT;
		open OUT,">$outdir/$t[1].eventalign.txt";
	}
	print OUT "$t[0]\t$t[2]\t$t[3]\t$t[4]\n";
	$kmer=$t[1];
}
close IN; close OUT;
