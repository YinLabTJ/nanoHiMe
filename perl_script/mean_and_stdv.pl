$input_file=shift; $outdir=shift;
open IN,"$input_file";
my @t=split /\//,$input_file;
my $prefix=$t[-1];
`mkdir -p $outdir/$prefix\_mean_stdv`;
while(<IN>){
	chomp;
	@t=split;
	next if($t[9] eq "NNNNNN");
	if($read ne $t[3]){
		if(defined($read)){
			foreach $key(sort keys %sixmer){
				next if($count{$key}<=2);
				$mean=$total{$key}/$count{$key};
				$n2=0;
				@sample=split /,/,$sixmer{$key};
				foreach $sample(@sample){
					$n2+=(($sample-$mean)*($sample-$mean));
				}
				$stdv=$n2/($count{$key}-1);
				$stdv=$stdv**0.5;
				$mean=sprintf "%.2f",$mean;
				$stdv=sprintf "%.2f",$stdv;
				print OUT "$key\t$mean\t$stdv\t$count{$key}\t$model_mean{$key}\t$model_stdv{$key}\n";
			}
		}
		%sixmer=();  %total=(); %count=();
		$read=$t[3];
		close OUT;
		open OUT,">$outdir/$prefix\_mean_stdv/$read.mean_and_stdv_out";
	}
	$sixmer{$t[9]}.="," if($sixmer{$t[9]});
	$sixmer{$t[9]}.=$t[6];
	$total{$t[9]}+=$t[6];
	$count{$t[9]}++;
	$model_mean{$t[9]}=$t[10];
	$model_stdv{$t[9]}=$t[11];
}
close IN;
foreach $key(sort keys %sixmer){
	next if($count{$key}<=2);
	$mean=$total{$key}/$count{$key};
	$n2=0;
	@sample=split /,/,$sixmer{$key};
	foreach $sample(@sample){
		$n2+=(($sample-$mean)*($sample-$mean));
	}
	$stdv=$n2/($count{$key}-1);
	$stdv=$stdv**0.5;
	$mean=sprintf "%.2f",$mean;
	$stdv=sprintf "%.2f",$stdv;
	print OUT "$key\t$mean\t$stdv\t$count{$key}\t$model_mean{$key}\t$model_stdv{$key}\n";
}
close OUT;
