#conda activate circos
#cp ../Can.genome ./Can_karyotype.txt
#cp ../Col.genome ./Col_karyotype.txt
# manually edit

#minimap2 -cx asm5 -t 10 ../Can_final_23022022_new_hapog.fasta ../Col_final_23022022_new_hapog_patch_hapog_reorder.fasta > Col_Can_align.paf

#awk '{if ($5 == "-") print $1,$3,$4,$6,$9,$8,$12,"identity="$10/$11; else print $1,$3,$4,$6,$8,$9,$12,"identity="$10/$11}' Col_Can_align.paf > Col_Can.links # last column is the identity

awk '$7>30' Col_Can.links | awk '{print $1,$2,$3,$4,$5,$6,$8}'> Col_Can_unique.links

awk '{print $1,$2,$3,$4,$5,$6,$8}' Col_Can.links > Col_Can_non-uniq.links

/tmp/Ziming_analysis/synteny_Col_Can/circos/circos-tools-0.23/tools/bundlelinks/bin/bundlelinks -links Col_Can_unique.links -max_gap 1e4 > Col_Can_unique.links.bundles

/tmp/Ziming_analysis/synteny_Col_Can/circos/circos-tools-0.23/tools/bundlelinks/bin/bundlelinks -links Col_Can_non-uniq.links -max_gap 1e4 > Col_Can_non-uniq.links.bundles

awk '$5 == "-"' Col_Can_align.paf | awk '{print $1,$3,$4,$6,$9,$8,$12,"identity="$10/$11}' | awk '$7>30' | awk '{print $1,$2,$3,$4,$5,$6}' > Col_Can_inverted.link

awk '$5 == "-"' Col_Can_align.paf | awk '{print $1,$3,$4,$6,$9,$8,$12,"identity="$10/$11}' | awk '{print $1,$2,$3,$4,$5,$6}' > Col_Can_non-uniq_inverted.link

/tmp/Ziming_analysis/synteny_Col_Can/circos/circos-tools-0.23/tools/bundlelinks/bin/bundlelinks -links Col_Can_inverted.link -max_gap 1e4 > Col_Can_inverted.links.bundles

/tmp/Ziming_analysis/synteny_Col_Can/circos/circos-tools-0.23/tools/bundlelinks/bin/bundlelinks -links Col_Can_non-uniq_inverted.link -max_gap 1e4 > Col_Can_non-uniq_inverted.links.bundles


# need to set the identity to a range of [1-6]
#awk '$7>30' Col_Can.links | awk '{if ($8 > 0.95) {$9="1"} else if ($8 > 0.90) {$9="2"} else if ($8 > 0.80) {$9="3"} else {$9="4"} print;}' | awk '$9="value="$9' | awk '{print $1,$2,$3,$4,$5,$6,$9}'> Col_Can_unique_good.links

#minimap2 -cx asm5 -t 10 ../Col_final_23022022_new_hapog_patch_hapog_reorder.fasta ../Col_final_23022022_new_hapog_patch_hapog_reorder.fasta > Col_col_align.paf
#awk '{print $1,$3,$4,$6,$8,$9,$5,$10,$11,$12}' Col_col_align.paf | awk '$11=$8/$9' > Col_Col.links
#minimap2 -cx asm5 -t 10 ../Can_final_23022022_new_hapog.fasta ../Can_final_23022022_new_hapog.fasta > Can_can_align.paf
#awk '{print $1,$3,$4,$6,$8,$9,$5,$10,$11,$12}' Can_can_align.paf | awk '$11=$8/$9' > Can_Can.links

# organized the Telomere repeat motif search file locally
#sed -i 's/Chr/coChr/g' Col_telomeric_repeat_windows.txt
#sed -i 's/Chr/caChr/g' Can_telomeric_repeat_windows.txt 

#awk -F "\t" '($4>10) && ($5>10)' Can_telomeric_repeat_windows.txt | sed 's/\t/ /g' | awk '{print $1,$2,$3,$4","$5","$6}' > Can_telomeric_selected.txt

#awk -F "\t" '($4>10) && ($5>10)' Col_telomeric_repeat_windows.txt | sed 's/\t/ /g' | awk '{print $1,$2,$3,$4","$5","$6}' > Col_telomeric_selected.txt

#sed 's/\t/ /g' Col_telomeric_repeat_windows.txt | awk '{print $1,$2,$3,$4","$5","$6}' > Col_telomeric.txt
#sed 's/\t/ /g' Can_telomeric_repeat_windows.txt | awk '{print $1,$2,$3,$4","$5","$6}' > Can_telomeric.txt

#cat Col_telomeric.txt Can_telomeric.txt > Col_Can_telomeric_repeat.txt
#awk '{print $1,$2,$3,$4","$5","$6}' Col_Can_telomeric_repeat.txt > Col_Can_telomeric_repeat1.txt
#mv Col_Can_telomeric_repeat1.txt Col_Can_telomeric_repeat.txt

#awk -F "\t" '{print $1,$2,$3,$11}' Col.modified_bases.5mC.bed | sed 's/Chr/coChr/g' > Col.CG.bed
#awk -F "\t" '{print $1,$2,$3,$11}' Can.modified_bases.5mC.bed | sed 's/Chr/caChr/g' > Can.CG.bed

#bedtools makewindows -g ../Col.genome -w 10000 > Col_1k.bed
#bedtools makewindows -g ../Can.genome -w 10000 > Can_1k.bed
#sed -i 's/Chr/caChr/g' Can_1k.bed
#sed -i 's/Chr/coChr/g' Col_1k.bed
#awk '$2=$2+1' Col_1k.bed | sort -k1,1 -k2,2n > Col_1k_new.bed
#awk '$2=$2+1' Can_1k.bed | sort -k1,1 -k2,2n > Can_1k_new.bed
#sed -i 's/ /\t/g' Col_1k_new.bed
#sed -i 's/ /\t/g' Can_1k_new.bed
#sort -k1,1 -k2,2n Col.CG.bed > Col.CG.sorted.bed
#sort -k1,1 -k2,2n Can.CG.bed > Can.CG.sorted.bed
#sed -i 's/ /\t/g' Can.CG.sorted.bed
#bedtools map -a Col_1k_new.bed -b Col.CG.sorted.bed -c 4 -o mean > Col.CG.1k.bed
#bedtools map -a Can_1k_new.bed -b Can.CG.sorted.bed -c 4 -o mean > Can.CG.1k.bed

#cat Col_telomeric_selected.txt Can_telomeric_selected.txt  > Can_Col_telomeric_selected.txt
#cat Col.CG.1k.bed Can.CG.1k.bed > Can_Col.CG.1k.bed

sort -k1,1 -k2,2n /tmp/Ziming_analysis/synteny_Col_Can/genomecov/Col.hifi.coverage > Col.hifi.coverage.sorted
sort -k1,1 -k2,2n /tmp/Ziming_analysis/synteny_Col_Can/genomecov/Can.hifi.coverage > Can.hifi.coverage.sorted
#sort -k1,1 -k2,2n /tmp/Ziming_analysis/synteny_Col_Can/genomecov/Col.ont.coverage > Col.ont.coverage.sorted
#sort -k1,1 -k2,2n /tmp/Ziming_analysis/synteny_Col_Can/genomecov/Can.ont.coverage > Can.ont.coverage.sorted

bedtools map -a Col_1k_new.bed -b Col.hifi.coverage.sorted -c 4 -o mean > Col.hifi.coverage.1k.bed

bedtools map -a Can_1k_new.bed -b Can.hifi.coverage.sorted -c 4 -o mean > Can.hifi.coverage.1k.bed

#bedtools map -a Col_1k_new.bed -b Col.ont.coverage.sorted -c 4 -o mean > Col.ont.coverage.1k.bed

#bedtools map -a Can_1k_new.bed -b Can.ont.coverage.sorted -c 4 -o mean > Can.ont.coverage.1k.bed

cat Col.hifi.coverage.1k.bed Can.hifi.coverage.1k.bed > Can_Col.hifi.coverage.1k.bed
cat Col.ont.coverage.1k.bed Can.ont.coverage.1k.bed > Can_Col.ont.coverage.1k.bed

#conda activate circos
#circos --conf circos.conf

/home/ucbtzz1/seqtk/seqtk cutN -n 5 -g ../Can_final_23022022_new_hapog.fasta > Can_gap.txt
/home/ucbtzz1/seqtk/seqtk cutN -n 5 -g ../Col_final_23022022_new_hapog_patch_hapog_reorder.fasta > Col.gap.txt
cat Col.gap.txt Can_gap.txt > Can_Col.gap.txt


cp ../intervals/Can_GCcontent_10k.bed ./
cp ../intervals/Col_GCcontent_10k.bed ./

sed -i 's/Chr/caChr/g' Can_GCcontent_10k.bed
sed -i 's/Chr/coChr/g' Col_GCcontent_10k.bed 

cat Col_GCcontent_10k.bed Can_GCcontent_10k.bed > Col_Can_GCcontent_10k.bed

cat ../Col_repeat_density.txt ../Can_repeat_density.txt > Col_Can_repeat_density


