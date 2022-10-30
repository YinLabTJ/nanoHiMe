use FindBin;
use File::Basename;

my $perl_main="$FindBin::Bin/$FindBin::Script";
my ($filename, $program_dir) = fileparse($perl_main);

$kmer=shift;
$kmer_dir=shift;
$filter_file=shift;
$abd_file=shift;
$iteration_limit=shift;
$delta_mu1=shift;
$delta_sigma1=shift;
$omega1=shift;
$delta_mu2=shift;
$delta_sigma2=shift;
$omega2=shift;


$e=2.718281828459;
$pi=3.1415926535898;

$delta_mu1 ||= 0;
$delta_sigma1 ||=0;
$delta_mu2 ||= 0;
$delta_sigma2 ||=1;

$omega1 ||= 0.2;
$omega2 ||= 0.8;

unless($omega1+$omega2){
	$omega1 = 0.2; $omega2 = 0.8; 
}

open IN,"$filter_file";
while(<IN>){
	chomp;
	@t=split /\//;
	$t[-1]=~s/\.mean_and_stdv_out//g;
	$num=<IN>; chomp($num);
	$filter{$t[-1]}=1 if($num<30);
}
close IN;

open IN,"$program_dir/../model/PCR/parameter.sort.txt";
while(<IN>){
	chomp;
	@t=split;
	if($t[0] eq $kmer){
		$model_mu=$t[1];
		$model_sigma=$t[2];
	}
}
close IN;

open IN,"$abd_file";
while(<IN>){
	chomp;
	@t=split;
	$a{$t[0]}=$t[1];
	$b{$t[0]}=$t[2];
	$d{$t[0]}=$t[3];
}
close IN;

open IN,"$kmer_dir/$kmer.eventalign.txt";
$i=0;
while(<IN>){
	chomp;
	@t=split;
	next if($filter{$t[0]});
	$E[$i]=$t[1];
	if($b{$t[0]}>0){
		$para_sigma[$i]=$d{$t[0]}/$b{$t[0]};
	}else{
		$para_sigma[$i]=$d{$t[0]}/$b{$t[0]};
	}
	$i++;	
}
close IN;


#methylated initial value
$mu1=$model_mu+$delta_mu1; $sigma1=$model_sigma+$delta_sigma1;
$mu2=$model_mu+$delta_mu2; $sigma2=$model_sigma+$delta_sigma2;


my $round=0;
while(1){
	#step E
	last if(@E<50);
	for($i=0;$i<@E;$i++){
		$P1=$omega1/((2*$pi*($sigma1*$para_sigma[$i])**2)**0.5)*$e**(-1*(($E[$i]-$mu1)**2)/(2*($sigma1*$para_sigma[$i])**2));
		$P2=$omega2/((2*$pi*($sigma2*$para_sigma[$i])**2)**0.5)*$e**(-1*(($E[$i]-$mu2)**2)/(2*($sigma2*$para_sigma[$i])**2));
		$Rj1[$i]=$P1/($P1+$P2);
		$Rj2[$i]=$P2/($P1+$P2);
		$Pt=$Pt+$P1+$P2;
	}
	#step M
	$sumRj1=0; $sumRj2=0; $sumRj1e=0; $sumRj2e=0; $sumRj1sigma=0; $sumRj2sigma=0;  
	for($i=0;$i<@E;$i++){
		$sumRj1+=$Rj1[$i];
		$sumRj2+=$Rj2[$i];
		$sumRj1e+=$Rj1[$i]*$E[$i];
		$sumRj2e+=$Rj2[$i]*$E[$i];
	}
	$omega1=$sumRj1/@E;
	$omega2=$sumRj2/@E;
	$mu1=$sumRj1e/$sumRj1;
	$mu2=$sumRj2e/$sumRj2;
	for($i=0;$i<@E;$i++){
		$sumRj1sigma+=$Rj1[$i]*(($E[$i]-$mu1)/$para_sigma[$i])**2;
		$sumRj2sigma+=$Rj2[$i]*(($E[$i]-$mu2)/$para_sigma[$i])**2;
	}
	$sigma1=($sumRj1sigma/$sumRj1)**0.5;
	$sigma2=($sumRj2sigma/$sumRj2)**0.5;

	
	if(abs($lastmu1-$mu1)<0.0001 && abs($lastmu2-$mu2)<0.0001){
		print "$kmer\t$omega1\t$mu1\t$sigma1\t$omega2\t$mu2\t$sigma2\n";
		last;
	}
	if($sigma1<0.00001 || $sigma2<0.00001){
		print "$kmer\t$omega1\t$mu1\t$sigma1\t$omega2\t$mu2\t$sigma2\n";
		last;
	}
	$round++;
	if(defined($iteration_limit) && $iteration_limit<=$round){
		print "$kmer\t$omega1\t$mu1\t$sigma1\t$omega2\t$mu2\t$sigma2\n";
		last;
	}
	$lastmu1=$mu1; $lastmu2=$mu2;
}
