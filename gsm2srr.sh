#!/bin/bash
# entrez-direct  （运行需要依赖这个包） conda install entrez-derect
#gse98149.k9.(gsm2sample.txt) 输入文件后缀
ls GSM_to_SRR_encode.txt &&  rm GSM_to_SRR_encode.txt


GSMs=`cat $1| cut -f1 `
input=$1
for GSM in $GSMs
do
 SRR_string=`esearch -db sra -query $GSM | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc`
 #echo -e "$GSM\t$SRR_string" >> GSM_to_SRR_tsv.txt
 SRR_array_unsorted=($SRR_string) #将$SRR_string 放入一对括号内，就变成数组
 IFS=$'\n' SRR_array=($(sort <<<"${SRR_array_unsorted[*]}"))  #${SRR_array_unsorted[*]} 将数组里所有元素进行排序，放到一堆括号里，所以排序后的也是一个数组，命名为SRR_array
 unset IFS
 
 #srr_total=${#SRR_array[@]} #的意思是获取数组SRR_array[@] 的个数
 for SRR in ${SRR_array[@]}
 do
 echo -e "$GSM\t$SRR" >> GSM_to_SRR_encode.txt
 done
done


awk -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{print $2,a[$1]}' $input GSM_to_SRR_encode.txt >${input/gsm2sample/srr2sample}
