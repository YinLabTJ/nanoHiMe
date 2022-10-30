for  i in `find $1 -name "*.mean_and_stdv_out"`
do
	j=${i/\.mean_and_stdv_out/}
	mkdir -p $2
	grep -v CG $i | grep -v ^.....C > $2/coefficients.input
	cd $2/
	Rscript $3/Rscript/coefficients.R  > coefficients.out
	perl $3/perl_script/coefficients.pl $j >> $2/coefficients.abd
done
