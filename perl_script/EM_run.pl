use FindBin;
use File::Basename;

my $perl_main="$FindBin::Bin/$FindBin::Script";
my ($filename, $dir) = fileparse($perl_main);

$kmer_parameter_file=shift;
$kmer_dir=shift;
$filter_file=shift;
$abd_file=shift;
$gaussian_num=shift;
$iteration_EM=shift;
$EM_shift_file=shift;
$output_file=shift;

$gaussian_num ||= 3;
$iteration_EM ||= 1000;
$EM_shift_file ||= "$dir/EM_shift.txt";
$output_file ||= "./EM.out";


open IN,"$EM_shift_file";
<IN>;
while(<IN>){
	chomp;
	push @EM_para,$_;
}
close IN;

if(-e $output_file){
	`mv $output_file $output_file.backup`;
}
open OUT,">$output_file";
print OUT "K-mer";
for($i=1;$i<=$gaussian_num;$i++){
	print OUT "\tomega_$i\tmu_$i\tsigma_$i";
}
print OUT "\n";
close OUT;
open IN,"$kmer_parameter_file";
while(<IN>){
	chomp;
	@t=split;
	foreach $EM_p(@EM_para){
		`perl $dir/EM.$gaussian_num.pl $t[0] $kmer_dir $filter_file $abd_file $iteration_EM $EM_p >> $output_file`;
	}
}
close IN; close SH;
