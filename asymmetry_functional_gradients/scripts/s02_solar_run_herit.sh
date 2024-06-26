#! /bin/sh
# this is a universal ICBM processing script

if [ $# -ne 2 ]
then
        echo $1
        echo Give the name for the genetics header and data files  such as out.header and out.txt
        exit
fi



# don't modify after this line

pheno_list=`more $1`
echo  $

# Lets make solar_run file
out_file=$2'_'solar.run.inorm
echo pheno load $2 >  $out_file
echo model new >>$out_file
echo covar age^1,2#sex >>$out_file


for cur_dir  in  $pheno_list
do
echo define $cur_dir'_'INORM = inorm_$cur_dir >> $out_file

echo trait $cur_dir'_'INORM >> $out_file
echo polyg -all -s >>$out_file

done
