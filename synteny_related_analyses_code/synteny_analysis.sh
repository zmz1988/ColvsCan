#conda activate omni-C

#tidk find -f /tmp/Ziming_analysis/annotation/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta -c plants -o Col
#tidk find -f /tmp/Ziming_analysis/annotation/Can_final_23022022_new_hapog.fasta -c plants -o Can
#tidk plot -c finder/Can_telomeric_repeat_windows.csv -o Can
#tidk plot -c finder/Col_telomeric_repeat_windows.csv -o Col

# build all kinds of intervals (https://wiki.bits.vib.be/index.php/Create_a_GC_content_track)
#cut -f 1,2 ../annotation/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta.fai > Col.genome
#cut -f 1,2 ../annotation/Can_final_23022022_new_hapog.fasta.fai > Can.genome
width=1000
#bedtools makewindows -g Col.genome -w ${width} > Col_${width}bps.bed
#bedtools makewindows -g Can.genome -w ${width} > Can_${width}bps.bed 
#bedtools nuc -fi ../annotation/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta -bed Col_${width}bps.bed > Col_nuc_${width}bps.txt
#bedtools nuc -fi ../annotation/Can_final_23022022_new_hapog.fasta -bed Can_${width}bps.bed > Can_nuc_${width}bps.txt
# width was defined earlier and was set to 1000, hence the title column filled with 'GCpc_1000bps'
#gawk -v w=${width} 'BEGIN{FS="\t"; OFS="\t"} {if (FNR>1) {print $1,$2,$3,"GCpc_"w"bps",$5}}' Col_nuc_${width}bps.txt > Col_${width}bps.igv
#gawk -v w=${width} 'BEGIN{FS="\t"; OFS="\t"} {if (FNR>1) {print $1,$2,$3,"GCpc_"w"bps",$5}}' Can_nuc_${width}bps.txt > Can_${width}bps.igv    
#igvtools toTDF -z 5 -f min,max,mean Col_${width}bps.igv Col_${width}bps.tdf ../annotation/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta
#igvtools toTDF -z 5 -f min,max,mean Can_${width}bps.igv Can_${width}bps.tdf ../annotation/Can_final_23022022_new_hapog.fasta

cut -f 1,2,3,5 Col_nuc_${width}bps.txt > Col_GC.bed
cut -f 1,2,3,5 Can_nuc_${width}bps.txt > Can_GC.bed


width=10000
bedtools makewindows -g Col.genome -w ${width} > Col_${width}bps.bed
bedtools makewindows -g Can.genome -w ${width} > Can_${width}bps.bed 
bedtools nuc -fi ../annotation/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta -bed Col_${width}bps.bed > Col_nuc_${width}bps.txt
bedtools nuc -fi ../annotation/Can_final_23022022_new_hapog.fasta -bed Can_${width}bps.bed > Can_nuc_${width}bps.txt

cut -f 1,2,3,5 Col_nuc_${width}bps.txt > Col_GCcontent_10k.bed
cut -f 1,2,3,5 Can_nuc_${width}bps.txt > Can_GCcontent_10k.bed

grep -v "Parent" /tmp/Ziming_analysis/annotation/EDTA_RepbaseLib/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta.mod.EDTA.TEanno.gff3 | sed 's/_p/_polished/g' > Col_final_23022022_new_hapog_patch_hapog_reorder.fasta.mod.EDTA.TEanno_onlyParent.gff3

grep -v "Parent" /tmp/Ziming_analysis/annotation/EDTA_RepbaseLib/Can_final_23022022_new_hapog.fasta.mod.EDTA.TEanno.gff3 | sed 's/_p/_polished/g' > Can_final_23022022_new_hapog.fasta.mod.EDTA.TEanno_onlyParent.gff3

bedtools intersect -a Col_${width}bps.bed -b Col_final_23022022_new_hapog_patch_hapog_reorder.fasta.mod.EDTA.TEanno_onlyParent.gff3 -c > Col_repeat_density.txt
bedtools intersect -a Can_${width}bps.bed -b Can_final_23022022_new_hapog.fasta.mod.EDTA.TEanno_onlyParent.gff3 -c > Can_repeat_density.txt
sed -i 's/Chr/caChr/g' Can_repeat_density.txt
sed -i 's/Chr/coChr/g' Col_repeat_density.txt


