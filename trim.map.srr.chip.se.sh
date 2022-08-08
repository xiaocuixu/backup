#!/bin/bash

bt2_index=~/01.ref/mm9Bowtie2Index/mm9
fastq_dir=download
cat sampleinfo/gse98149/gse98149.k9.srr2sample.txt  |sort -k1,1 |datamash -g 2 collapse 1 | while read line
do
units=(`echo $line`)
srr_str=${units[1]}
sample=${units[0]}

#:<<!
fn_arr=() #定义一个空数组，用来存放不同的srr trim之后的文件名
oldIFS=$IFS
IFS=,
srr_arr=($srr_str) #加括号变成数组
	for ((i=0;i<${#srr_arr[@]};i++))
	do
	#如果是双端数据则重命名
	ls $fastq_dir/${srr_arr[$i]}_2.fastq.gz && mv $fastq_dir/${srr_arr[$i]}_1.fastq.gz $fastq_dir/${srr_arr[$i]}.fastq.gz
	#对每一个srr做trim
	trim_galore -j 7 --basename ${srr_arr[$i]}  $fastq_dir/${srr_arr[$i]}.fastq.gz 
	fn_arr[$i]="${srr_arr[$i]}_trimmed.fq.gz" #trim之后的文件名
	#删除report
	rm ${srr_arr[$i]}.fastq.gz_trimming_report.txt
	done
IFS=$oldIFS


fn=$(IFS=,;echo "${fn_arr[*]}") #将不同srr连成一个字符串
bowtie2 -p 30 -x $bt2_index -U $fn -S $sample.sam
samtools sort -@ 30 -o bam/$sample.bam $sample.sam
rm ${fn_arr[@]} $sample.sam

samtools view -F 4 -o bam/$sample.F4.bam bam/$sample.bam
java -jar ~/02.bin/picard.jar MarkDuplicates \
      I=bam/$sample.F4.bam \
      O=bam/$sample.DeDup.bam \
      REMOVE_DUPLICATES=true \
      M=$sample.marked_dup_metrics.txt

rm $sample.marked_dup_metrics.txt 
samtools flagstat bam/$sample.bam >flagstat/$sample.flagstat
samtools flagstat bam/$sample.DeDup.bam >flagstat/$sample.DeDup.flagstat
rm bam/$sample.F4.bam #bam/$sample.bam
samtools index bam/$sample.DeDup.bam
bamCoverage -p max/2 --bam bam/$sample.DeDup.bam -o bw/$sample.RPKM.bw --normalizeUsing RPKM

all=`awk 'NR==1{print $1}' flagstat/${sample}.flagstat`
map=`awk 'NR==5{print $1}' flagstat/${sample}.flagstat`
dedup=`awk 'NR==5{print $1}' flagstat/${sample}.DeDup.flagstat`
awk 'BEGIN{OFS="\t";print "'$sample'","'$all'","'$map'","'$map'"/"'$all'","'$dedup'","'$dedup'"/"'$all'"}' >>SEReadsCount.txt
#!
rm bam/$sample.bam
done


