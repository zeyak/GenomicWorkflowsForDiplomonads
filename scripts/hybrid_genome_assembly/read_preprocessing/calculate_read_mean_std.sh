#MaSuRCA requires mean and stdev for the reads

#Calculation of the 'mean and stdev' of Reads: (-b binary parameter)
#Run for each pair reads: run1, run2, run3
for i in path_to_the_pair_reads/run1*
do
echo  "$i   "
awk -b 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sqrt(sq/n-m*m));}' $i
done