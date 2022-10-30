for  i in `find $1/ -name "*.mean_and_stdv_out"`
do
	echo "$i" >> $1/filtered.num.txt
	grep -v A $i | wc -l >> $1/filtered.num.txt
done
