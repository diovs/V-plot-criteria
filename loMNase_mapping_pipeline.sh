
ls *_R1.fq.gz |while read id;
do
filename=`basename $id _R1.fq.gz`
#trim adapter
raw=`zcat ${id} |wc -l`
read_count=`echo "${raw}/4" | bc`
echo "raw: ${read_count} (100%)" | tee -a ${filename}_Report.txt

trim_galore -j 80 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired ${filename}_R1.fq.gz ${filename}_R2.fq.gz --gzip
#fastp -w 16 -i ${filename}_R1_val_1.fq.gz -I ${filename}_R2_val_2.fq.gz -o ${filename}_R1.QC.fq.gz -O ${filename}_R2.QC.fq.gz
###reads mapping
trim=`zcat ${filename}_R1_val_1.fq.gz |wc -l`
trim_read_count=`echo "${trim}/4" | bc`
QC_rate=`echo "scale=2; 100*${trim_read_count}/${read_count}" | bc`
echo "QC: ${trim_read_count} (${QC_rate}%)" | tee -a ${filename}_Report.txt

fastqc -t 80 --nogroup ${filename}_R?_val_?.fq.gz
mv fastp.html ${filename}'_'fastp.html
mv *html fastq
rm *json *zip
bowtie2 -t -p 60 --no-unal --no-mixed --no-discordant --dovetail --very-sensitive --score-min L,0,-0.4 -X 1000 \
-x /ssd/index/bowtie2/hg38XY+mm10XY -1 ${filename}_R1_val_1.fq.gz -2 ${filename}_R2_val_2.fq.gz | samtools view -bh -q 10 -@ 80 -  > ${filename}'_'q10.bam


bamToBed -bedpe -mate1 -i ${filename}_q10.bam > ${filename}_q10.bed

mapping_q10_count=`cat ${filename}_q10.bed |wc -l`
mapping_q10_rate=`echo "scale=2; 100*${mapping_q10_count}/${trim_read_count}" | bc`
echo "mapping & q10: ${mapping_q10_count} (${mapping_q10_rate}%)" | tee -a ${filename}_Report.txt

awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' ${filename}_q10.bed | sort -k1,1 -k2,2n -k3,3n -k4,4 -u -S 100% | cut -f 1-3 > ${filename}_q10_rmdup.bed
dedup_count=`cat ${filename}_q10_rmdup.bed |wc -l`
dedup_rate=`echo "scale=2; 100*${dedup_count}/${mapping_q10_count}" | bc`
echo "dedup: ${dedup_count} (${dedup_rate}%)" | tee -a ${filename}_Report.txt

awk '$1!~/chr[CLMT]/' OFS='\t' ${filename}_q10_rmdup.bed > ${filename}_q10_rmdup_rmchrM.bed
rm_chrM_count=`cat ${filename}_q10_rmdup_rmchrM.bed |wc -l`
rm_chrM_rate=`echo "scale=2; 100*${rm_chrM_count}/${dedup_count}" | bc`
echo "rm_chrM: ${rm_chrM_count} (${rm_chrM_rate}%)" | tee -a ${filename}_Report.txt


#awk '{print $3-$2}' OFS='\t' ${filename}.bed > ${filename}_frag_length.bed


awk '{if($3-$2<=100){print $0}}' OFS='\t' ${filename}_q10_rmdup_rmchrM.bed > ${filename}_q10_rmdup_rmchrM_frag_filter.bed
target_fragment_count=`cat ${filename}_q10_rmdup_rmchrM_frag_filter.bed|wc -l`

target_fragment_rate=`echo "scale=2; 100*${target_fragment_count}/${rm_chrM_count}" | bc`
echo "target_fragment_rate: ${target_fragment_count} (${target_fragment_rate}%)" | tee -a ${filename}_Report.txt

awk '$1!~/^Mm/' OFS='\t' ${filename}_q10_rmdup_rmchrM_frag_filter.bed > ${filename}_q10_rmdup_rmchrM_frag_filter_rm3T3.bed
awk '$1 ~/^Mm/' OFS='\t' ${filename}_q10_rmdup_rmchrM_frag_filter.bed  > ${filename}_3T3.bed


NIH3T3_read_count=`cat ${filename}_3T3.bed|wc -l`
NIH3T3_rate=`echo "scale=2; 100*${NIH3T3_read_count}/${target_fragment_count}" | bc`
echo "3T3_rate: ${NIH3T3_read_count} (${NIH3T3_rate}%)" | tee -a ${filename}_Report.txt

rm3T3_read_count=`cat ${filename}_q10_rmdup_rmchrM_frag_filter_rm3T3.bed | wc -l`

genomeCoverageBed -bg -i ${filename}_q10_rmdup_rmchrM_frag_filter_rm3T3.bed -g /ssd/genome/hg38_chromsize.txt | awk -v a=${rm3T3_read_count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | awk '{$4/=1;print}' OFS='\t' | sort -k1,1 -k2,2n -S 100% > ${filename}.bdg
bedGraphToBigWig ${filename}.bdg /ssd/genome/hg38_chromsize.txt ${filename}.bw


sed s/Mm.// ${filename}_3T3.bed|genomeCoverageBed -bg -i - -g /ssd/genome/mm10_chromsize.txt | awk -v a=${NIH3T3_read_count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | awk '{$4/=1;print}' OFS='\t' | sort -k1,1 -k2,2n -S 100% > ${filename}_Mm.bdg
bedGraphToBigWig ${filename}_Mm.bdg /ssd/genome/mm10_chromsize.txt ${filename}_Mm.bw


rm ${filename}_R1_val_1.fq.gz ${filename}_R2_val_2.fq.gz ${filename}_q10.bed ${filename}_q10_rmdup.bed \
 ${filename}_q10_rmdup_rmchrM.bed ${filename}_q10_rmdup_rmchrM_frag_filter.bed
rm *.bdg
mkdir ${filename}
mv ${filename}* ${filename}/
done

