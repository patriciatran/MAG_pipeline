cp results/SRR9162907/checkm/quality_summary.tsv results/SRR9162907/checkm_for_dRep/quality_summary.csv 
sed -i 's|Bin Id|genome|g' results/SRR9162907/checkm_for_dRep/quality_summary.csv 
sed -i 's|Completeness|completeness|g' results/SRR9162907/checkm_for_dRep/quality_summary.csv  
sed -i 's|Contamination|contamination|g' results/SRR9162907/checkm_for_dRep/quality_summary.csv 
sed -i 's|\t|,|g' results/SRR9162907/checkm_for_dRep/quality_summary.csv 

awk -v OFS="," -F "," '$1=$1".fasta"' results/SRR9162907/checkm_for_dRep/quality_summary.csv > results/SRR9162907/checkm_for_dRep/quality_summary_awk.csv
sed -i 's|genome.fasta|genome|g' results/SRR9162907/checkm_for_dRep/quality_summary_awk.csv
#awk -v OFS="," -F "," '$1="results/SRR9162907/renamed_redined_MAGS/"$1' results/SRR9162907/checkm_for_dRep/quality_summary_awk.csv > results/SRR9162907/checkm_for_dRep/quality_summary_awk2.csv