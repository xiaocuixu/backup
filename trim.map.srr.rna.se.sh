#!/bin/bash

#公共数据的预处理，处理同一个样本有多个SRR的情况,单端mapping

#mkdir bw bam report fpkm_reps cufflinks_out
index=~/01.ref/mm9StarIndex/
fastq_dir=download
cat sampleinfo/gse76505/gse76505.rna.srr2sample.txt  |sort -k1,1 |datamash -g 2 collapse 1 |tail -n 15 |while read line
do
units=(`echo $line`)
srr_str=${units[1]}
sample=${units[0]}


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

fn=$(IFS=,;echo "${fn_arr[*]}") #将不同srr连成一个字符串,作为mapping软件的输入
STAR --runThreadN 30 --runMode alignReads --genomeDir $index --readFilesIn $fn --outFileNamePrefix ./${sample}_ --outSAMtype BAM SortedByCoordinate  --outSAMstrandField intronMotif --readFilesCommand zcat

rm ${fn_arr[@]} ${sample}_Log.progress.out ${sample}_Log.out
mv ${sample}_Aligned.sortedByCoord.out.bam bam

java -jar ~/02.bin/picard.jar MarkDuplicates \
      I=bam/${sample}_Aligned.sortedByCoord.out.bam \
      O=bam/$sample.DeDup.bam \
      REMOVE_DUPLICATES=true \
      M=$sample.marked_dup_metrics.txt

mv $sample.marked_dup_metrics.txt report
mv ${sample}_Log.final.out report
mv ${sample}_SJ.out.tab report

grep Unknown report/$sample.marked_dup_metrics.txt |awk '{print "'$sample'""\t"$10}' >QC.Dup.txt
grep 'Uniquely mapped reads number' report/${sample}_Log.final.out |awk '{print "'$sample'""\t"$6}'  >QC.Reads_num.txt
grep 'Uniquely mapped reads %' report/${sample}_Log.final.out |awk '{print "'$sample'""\t"$6}' >QC.Map_rate.txt
grep 'Number of input reads' report/${sample}_Log.final.out |awk '{print "'$sample'""\t"$6}' >QC.Input_Reads_num.txt
Rscript ~/03.script/qc.r
rm QC.Dup.txt QC.Reads_num.txt QC.Map_rate.txt QC.Input_Reads_num.txt

samtools index bam/${sample}_Aligned.sortedByCoord.out.bam
bamCoverage --bam bam/${sample}_Aligned.sortedByCoord.out.bam -o bw/$sample.RPKM.bw --normalizeUsing RPKM -p max/2

rm bam/$sample.DeDup.bam
cufflinks -G /home1/xuxc/ann/mm9.genes.chr.gtf -p 30 -o cufflinks_out/$sample bam/${sample}_Aligned.sortedByCoord.out.bam
awk 'NR>1{print $0}' cufflinks_out/${sample}/genes.fpkm_tracking |perl -ane 'if(! exists $h{$F[0]}){$h{$F[0]}=$_}else{$h{$F[0]}="NA\n"}END{foreach(keys %h){print $h{$_}}}' |grep -v '^NA' |cut -f1,10 |sort -k1b,1 >fpkm_reps/$sample.fpkm.uniq.tab

done


