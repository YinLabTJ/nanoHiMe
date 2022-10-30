$double=shift; $para=shift; $output=shift;
open IN,"$double";
<IN>;
while(<IN>){
	chomp;
	@t=split;
	$t[2]=sprintf "%.2f",$t[2];
	$t[3]=sprintf "%.2f",$t[3];
	$t[5]=sprintf "%.2f",$t[5];
	$t[6]=sprintf "%.2f",$t[6];
	$m1{$t[0]}=$t[2];  
	$std1{$t[0]}=$t[3];
	$m2{$t[0]}=$t[5];
	$std2{$t[0]}=$t[6];
}
close IN;

open IN,"$para";
open OUT,">$output";
while(<IN>){
	chomp;
	@t=split;
	if($m2{$t[0]}>0){
		print OUT "$t[0]\t$m1{$t[0]}\t$std1{$t[0]}\t$m2{$t[0]}\t$std2{$t[0]}\n";
	}elsif($m1{$t[0]}>0){
		print OUT "$t[0]\t$m1{$t[0]}\t$std1{$t[0]}\t$m1{$t[0]}\t$std1{$t[0]}\n";
	}else{
		print OUT "$t[0]\t$t[1]\t$t[2]\t$t[1]\t$t[2]\n";
	}
}
close IN; close OUT;
