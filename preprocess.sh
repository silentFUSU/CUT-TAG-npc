input=$1
peaks=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/data/Cut_tag/cleandata/D[0-6]/peaks/NPCD[0,1,2,3,6]_$1_merge/NPCD[0-6]_$1_peaks.narrowPeak)
cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript code/keep_regular_chroms.r  > $1/$1_merged_peaks.bed
bedtools intersect -a $1/$1_merged_peaks.bed -b  <(cat $peaks | awk -v OFS="\t" '{if ($10==-1){ print $1,int(($2+$3)/2),int(($2+$3)/2) }  else {print $1,$2+$10,$2+$10} }' ) -loj > $1/$1_merged_peaks.all_summits.txt
bash code/bed_to_saf.sh $1/$1_merged_peaks.bed $1/$1_merged_peaks.saf
echo "conversion to SAF done"
files=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/data/Cut_tag/cleandata/D[0-6]/bam/NPCD[0-6]_$1_[1-2].nodup.bam)
featureCounts -a $1/$1_merged_peaks.saf -o $1/$1_merged_peaks.counts $files -F SAF -T 8
