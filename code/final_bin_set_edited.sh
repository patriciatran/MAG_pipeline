awk -F "\t" '{ if(($12 >= 50) && ($13 <10)) {print} }' results/SRR9162907/checkm/quality_summary.tsv | cut -f 1 > results/SRR9162907/medium_quality_bins.txt
awk -F "\t" '{ if(($12 > 90) && ($13 <5)) {print} }' results/SRR9162907/checkm/quality_summary.tsv | cut -f 1 > results/SRR9162907/high_quality_bins.txt
cat results/SRR9162907/medium_quality_bins.txt results/SRR9162907/high_quality_bins.txt > results/SRR9162907/final_bin_set.txt
sort results/SRR9162907/final_bin_set.txt | uniq > results/SRR9162907/final_bin_set_unique.txt

mkdir results/SRR9162907/final_bin_set
sed -e "s|^|cp results/SRR9162907/dereplicated_bins/dereplicated_genomes/|g" results/SRR9162907/final_bin_set_unique.txt > results/SRR9162907/copy_final_bin_set.sh
sed -i "s|$|.fasta results/SRR9162907/final_bin_set/.|g" results/SRR9162907/copy_final_bin_set.sh 
bash results/SRR9162907/copy_final_bin_set.sh
