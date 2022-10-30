use FindBin;
use File::Basename;

my $perl_main="$FindBin::Bin/$FindBin::Script";
my ($filename, $dir) = fileparse($perl_main);


my $modification=shift; #ma, mcg or ma_mcg
my $round=shift;
my $ref_fasta=shift;
my $input_eventalign=shift;
my $outdir=shift;
my $gaussian_num=shift;
my $gaussian_params_file=shift;
my $EM_iteration=shift;

$round ||=5;
$gaussian_num ||=3;


$old_model="$dir/../model/PCR/parameter.sort.txt";
for($i=1;$i<=$round;$i++){
	`mkdir -p $outdir/round$i`;
	`$dir/../training/$modification/re_align.$modification $old_model $ref_fasta $input_eventalign > $outdir/round$i/training.eventalign.txt`;
	`perl $dir/mean_and_stdv.pl $outdir/round$i/training.eventalign.txt $outdir/round$i/`;
	`bash $dir/../shell_script/filter.$modification.sh $outdir/round$i/training.eventalign.txt_mean_stdv`;
	`bash $dir/../shell_script/abd.$modification.sh $outdir/round$i/training.eventalign.txt_mean_stdv $outdir/round$i/training.eventalign.txt_coefficients $dir/..`;
	`mkdir -p $outdir/round$i/6mer`;
	`$dir/../training/filter_border_check $outdir/round$i/training.eventalign.txt | sort -k2,2 > $outdir/round$i/training.eventalign.6mer.txt`;
	`perl $dir/split_kmer.pl $outdir/round$i/training.eventalign.6mer.txt $outdir/round$i/6mer`;
	`perl $dir/adjust.pl $outdir/round$i/training.eventalign.txt_coefficients/coefficients.abd $outdir/round$i/6mer $outdir/round$i/6mer_final`;
	`perl $dir/single_Gaussian.pl $outdir/round$i/training.eventalign.txt_mean_stdv $outdir/round$i/6mer_final $old_model | sort -k1,1 > $outdir/round$i/parameter.txt`;	     
	$old_model="$outdir/round$i/parameter.txt";
}

`perl $dir/EM_run.pl $old_model $outdir/round$round/6mer_final $outdir/round$round/training.eventalign.txt_mean_stdv/filtered.num.txt $outdir/round$round/training.eventalign.txt_coefficients/coefficients.abd $gaussian_num $EM_iteration $gaussian_params_file $outdir/EM.out`;
