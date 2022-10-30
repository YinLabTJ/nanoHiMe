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
$delta_mu3=shift;
$delta_sigma3=shift;
$omega3=shift;
$delta_mu4=shift;
$delta_sigma4=shift;
$omega4=shift;



$e=2.718281828459;
$pi=3.1415926535898;

$delta_mu1 ||= 0;
$delta_sigma1 ||=0;
$delta_mu2 ||= 0;
$delta_sigma2 ||=1;
$delta_mu3 ||= -5;
$delta_sigma3 ||=1;
$delta_mu3 ||= 5;
$delta_sigma3 ||=1;


$omega1 ||= 0.2;
$omega2 ||= 0.2;
$omega3 ||= 0.3;
$omega4 ||= 0.3;

unless($omega1+$omega2+$omega3+$omega4==1){
	$omega1 = 0.2; $omega2 = 0.2; $omega3 = 0.3; $omega4 = 0.3;
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
$mu3=$model_mu+$delta_mu3; $sigma3=$model_sigma+$delta_sigma3;
$mu4=$model_mu+$delta_mu4; $sigma4=$model_sigma+$delta_sigma4;


my $round=0;
while(1){
	#step E
	last if(@E<50);
	for($i=0;$i<@E;$i++){
		$P1=$omega1/((2*$pi*($sigma1*$para_sigma[$i])**2)**0.5)*$e**(-1*(($E[$i]-$mu1)**2)/(2*($sigma1*$para_sigma[$i])**2));
		$P2=$omega2/((2*$pi*($sigma2*$para_sigma[$i])**2)**0.5)*$e**(-1*(($E[$i]-$mu2)**2)/(2*($sigma2*$para_sigma[$i])**2));
		$P3=$omega3/((2*$pi*($sigma3*$para_sigma[$i])**2)**0.5)*$e**(-1*(($E[$i]-$mu3)**2)/(2*($sigma3*$para_sigma[$i])**2));
		$P4=$omega4/((2*$pi*($sigma4*$para_sigma[$i])**2)**0.5)*$e**(-1*(($E[$i]-$mu4)**2)/(2*($sigma4*$para_sigma[$i])**2));
		$Rj1[$i]=$P1/($P1+$P2+$P3+$P4);
		$Rj2[$i]=$P2/($P1+$P2+$P3+$P4);
		$Rj3[$i]=$P3/($P1+$P2+$P3+$P4);
		$Rj4[$i]=$P4/($P1+$P2+$P3+$P4)
		$Pt=$Pt+$P1+$P2+$P3+$P4;
	}
	#step M
	$sumRj1=0; $sumRj2=0; $sumRj1e=0; $sumRj2e=0; $sumRj1sigma=0; $sumRj2sigma=0;  
	$sumRj3=0; $sumRj3e=0; $sumRj3sigma=0; $sumRj4=0; $sumRj4e=0; $sumRj4sigma=0;
	for($i=0;$i<@E;$i++){
		$sumRj1+=$Rj1[$i];
		$sumRj2+=$Rj2[$i];
		$sumRj3+=$Rj3[$i];
		$sumRj4+=$Rj4[$i];
		$sumRj1e+=$Rj1[$i]*$E[$i];
		$sumRj2e+=$Rj2[$i]*$E[$i];
		$sumRj3e+=$Rj3[$i]*$E[$i];
		$sumRj4e+=$Rj4[$i]*$E[$i];
	}
	$omega1=$sumRj1/@E;
	$omega2=$sumRj2/@E;
	$omega3=$sumRj3/@E;
	$omega4=$sumRj4/@E;
	$mu1=$sumRj1e/$sumRj1;
	$mu2=$sumRj2e/$sumRj2;
	$mu3=$sumRj3e/$sumRj3;
	$mu4=$sumRj4e/$sumRj4;
	for($i=0;$i<@E;$i++){
		$sumRj1sigma+=$Rj1[$i]*(($E[$i]-$mu1)/$para_sigma[$i])**2;
		$sumRj2sigma+=$Rj2[$i]*(($E[$i]-$mu2)/$para_sigma[$i])**2;
		$sumRj3sigma+=$Rj3[$i]*(($E[$i]-$mu3)/$para_sigma[$i])**2;
		$sumRj4sigma+=$Rj4[$i]*(($E[$i]-$mu4)/$para_sigma[$i])**2;
	}
	$sigma1=($sumRj1sigma/$sumRj1)**0.5;
	$sigma2=($sumRj2sigma/$sumRj2)**0.5;
	$sigma3=($sumRj3sigma/$sumRj3)**0.5;
	$sigma4=($sumRj4sigma/$sumRj4)**0.5;
	
	if(abs($lastmu1-$mu1)<0.0001 && abs($lastmu2-$mu2)<0.0001 && abs($lastmu3-$mu3)<0.0001 && abs($lastmu4-$mu4)<0.0001){
		print "$kmer\t$omega1\t$mu1\t$sigma1\t$omega2\t$mu2\t$sigma2\t$omega3\t$mu3\t$sigma3\t$omega4\t$mu4\t$sigma4\n";
		last;
	}
	if($sigma1<0.00001 || $sigma2<0.00001 || $sigma3<0.00001 || $sigma4<0.00001){
		print "$kmer\t$omega1\t$mu1\t$sigma1\t$omega2\t$mu2\t$sigma2\t$omega3\t$mu3\t$sigma3\t$omega4\t$mu4\t$sigma4\n";
		last;	
	}
	$round++;
	if(defined($iteration_limit) && $iteration_limit<=$round){
		print "$kmer\t$omega1\t$mu1\t$sigma1\t$omega2\t$mu2\t$sigma2\t$omega3\t$mu3\t$sigma3\t$omega4\t$mu4\t$sigma4\n";
		last;
	}
	$lastmu1=$mu1; $lastmu2=$mu2; $lastmu3=$mu3; $lastmu4=$mu4;
}
