# (c) Patricia Tran
# About: This is a snakemake pipeline that will take your FASTQ reads and turn them into a set medium and high quality dereplicated MAGs, with taxonomy.

# Instructions:
# Change lines below to define input and output folders.
# Change individual settings for each step are in the "params" section of each "rule"

# Define input and output folders - Users, change this.
input_dir = "data/fastq"
output_dir = "results"

TMPDIR="/storage1/data10/tmp"

# from snakemake.io import glob_wildcards, for troubleshooting.
# Double brackets in {{sample}} for subfolders
fastq_files = glob_wildcards(f"{input_dir}/{{sample}}_R1.fastq.gz")
samples = list(fastq_files.sample)

## Rules below ##
rule all:
    input:
        # from sympy import expand
        expand("results/{sample}/taxonomy.tsv", sample=samples)
        
# Define rules for each step in the workflow
rule qc_reads:
    input:
        r1 = input_dir + "/{sample}_R1.fastq.gz",
        r2 = input_dir + "/{sample}_R2.fastq.gz"
    output:
        r1 = input_dir + "/qc/{sample}_R1.fastq.gz",
        r2 = input_dir + "/qc/{sample}_R2.fastq.gz",
        html = "results/{sample}/qc/{sample}_report.html",
        json = "results/{sample}/qc/{sample}_report.json"
    params:
        threads = 20
    conda: 
        "envs/fastp.yml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json} --thread {params.threads}
        """


rule assemble_reads:
# Description: This rules takes R1 and R2 paired FASTQ reads (QC-ed in the previous step), assuming they are metagenomic, and assembles them using SPADES.s
    input:
        r1 = input_dir + "/qc/{sample}_R1.fastq.gz",
        r2 = input_dir + "/qc/{sample}_R2.fastq.gz",
    output:
        assembly = directory(output_dir + "/{sample}/assembly")
    conda:
        "envs/spades.yml"
    params:
        memory = 100, #use max 600 on our server, unit is Gb
        kmers = "33,55,77,99,127",
        threads = 20

    shell:
        """
        spades.py -m {params.memory} -k {params.kmers} \
        --meta -t {params.threads} -1 {input.r1} -2 {input.r2} \
        -o {output.assembly}
        """

rule unzip_reads:
    input: 
        r1 = input_dir + "/qc/{sample}_R1.fastq.gz",
        r2 = input_dir + "/qc/{sample}_R2.fastq.gz",
        assembly = output_dir + "/{sample}/assembly"
    output:
        r1_unzip =  input_dir + "/unzipped/{sample}_1.fastq",
        r2_unzip =  input_dir + "/unzipped/{sample}_2.fastq"
    shell:
        """
        gzip -c -d {input.r1} > {output.r1_unzip}
        gzip -c -d {input.r1} > {output.r2_unzip}
        """

rule bin:
# Using the assembled scaffods (reads) in the assembly and the reads, we will now bin them using metawrap, using 3 different softwares.
    input:
        assembly = "results/{sample}/assembly/contigs.fasta",
        r1_unzip =  input_dir + "/unzipped/{sample}_1.fastq",
        r2_unzip =  input_dir + "/unzipped/{sample}_2.fastq"
    output:
        bins = directory("results/{sample}/bins/")
    conda:
        "envs/metawrap_env_2.yml"
    params:
        threads = 20,
        tmpdir = "storage1/data10/tmp"

    shell:
        """
        TMPDIR={params.tmpdir}
        metawrap binning -o {output.bins} -a {input.assembly} \
        -t {params.threads} \
        --metabat1 --metabat2 --maxbin2 \
        {input.r1_unzip} {input.r2_unzip}
        """

rule refine:
# Refine bins using dasTools. Requires to create a scaff to bin file first.
    input:
        bins = "results/{sample}/bins/",
        assembly = "results/{sample}/assembly/contigs.fasta",
        sample_dir = "results/{sample}/"
    output:
        refined_bins_tables = directory("results/{sample}/refine_bins_tables/"),
        refined_bins = directory("results/{sample}/refine_bins_DASTool_bins/")
    params:
        threads = 20
    conda:
        "envs/dasTool.yml"
    shell:
        """
        mkdir -p {output.refined_bins_tables}
        Fasta_to_Scaffolds2Bin.sh -e fa -i {input.bins}/maxbin2_bins/ > {output.refined_bins_tables}/maxbin2_scaf2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -e fa -i {input.bins}/metabat1_bins/ > {output.refined_bins_tables}/metabat1_scaf2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -e fa -i {input.bins}/metabat2_bins/ > {output.refined_bins_tables}/metabat2_scaf2bin.tsv
        
        DAS_Tool --write_bins 1 --write_bin_evals 1 --create_plots 1 \
        -i {output.refined_bins_tables}/metabat1_scaf2bin.tsv,{output.refined_bins_tables}/metabat2_scaf2bin.tsv,{output.refined_bins_tables}/maxbin2_scaf2bin.tsv \
        -l metabat1,metabat2,maxbin2 \
        -c {input.assembly} \
        -o {input.sample_dir}/refine_bins \
        --threads {params.threads}
        """

rule rename_MAGs:
# MAGs need to be called ".fasta" instead of ".fa" for dRep to work:
    input:
        refined_bins = "results/{sample}/refine_bins_DASTool_bins/",
        sample_dir = "results/{sample}"
    output:
        refined_bins_renamed = directory("results/{sample}/renamed_refined_MAGS/")
    shell:
        """
        mkdir -p {output.refined_bins_renamed}
        rename 's/.fa/.fasta/' {input.refined_bins}/*.fa
        mv {input.refined_bins}/*.fasta {output.refined_bins_renamed}/.      
        """

rule quality_check:
#Run checkm prior to running dRep because there is no way to set a tmp directory in dRep only. 
#We can then provide the .tsv file to dRep in the next step
     input:
         refined_bins_renamed = "results/{sample}/renamed_refined_MAGS/"
     output:
         folder = directory("results/{sample}/checkm/"),
         quality_summary = "results/{sample}/checkm/quality_summary.tsv"
     conda:
         "envs/checkM.yml"
     params:
         threads = 20,
         extension = "fasta", #file extension
         tmpdir = "/storage1/data10/tmp",
         pplacer_threads = 20
     shell:
        """
        checkm lineage_wf -x {params.extension} \
        -t {params.threads} \
        --tab_table \
        -f {output.quality_summary} \
        --tmpdir {params.tmpdir} \
        {input.refined_bins_renamed} {output.folder}
        """

rule edit_checkm_file:
#  # There are a few edits to the TSV file to be done --> convert to a CSV file, change the checkm headers Bin.Id to genome, Completeness to completeness, and Contamination to contamination
    input:
        quality_checkm = "results/{sample}/checkm/quality_summary.tsv" 
    output:
        checkm_for_drep = directory("results/{sample}/checkm_for_dRep/"),
        quality_summary_dRep = "results/{sample}/checkm_for_dRep/quality_summary.csv",
        quality_summary_dRep_awk = "results/{sample}/checkm_for_dRep/quality_summary_awk.csv"
    shell:
        """
        mkdir -p {output.checkm_for_drep}
        cp {input.quality_checkm} {output.quality_summary_dRep}
        sed -i 's|Bin Id|genome|g' {output.quality_summary_dRep}
        sed -i 's|Completeness|completeness|g' {output.quality_summary_dRep} 
        sed -i 's|Contamination|contamination|g' {output.quality_summary_dRep} 
        sed -i 's|\t|,|g' {output.quality_summary_dRep}

        # Finally, we need to add the file extension to the first MAG names in the first column:
        awk -v OFS="," -F "," '$1=$1".fasta"' {output.quality_summary_dRep} > {output.quality_summary_dRep_awk}
        sed -i 's|genome.fasta|genome|g' {output.quality_summary_dRep_awk}
        """

rule dereplicate:
# Once we have the bins, we will dereplicate them as to not have repetitive genomes. They need to be named ".fasta"
    input:
        refined_bins_renamed = "results/{sample}/renamed_refined_MAGS/",
        quality_summary_dRep_awk = "results/{sample}/checkm_for_dRep/quality_summary_awk.csv"
    output:
        dereplicated_bins = directory("results/{sample}/dereplicated_bins/"),
        bin_folder = "results/{sample}/dereplicated_bins/dereplicated_genomes/",
        quality_summary = "results/{sample}/dereplicated_bins/data_tables/genomeInfo.csv"
    conda:
        "envs/dRep.yml"
    params:
        completeness_min =  50, # Note that the defaut is 75
        sa = 0.99, # ANI threshold (default is 0.99)
        nc = 0.1, # coverage threshold (default is 0.1)
    shell:
        """
        dRep dereplicate \
        {output.dereplicated_bins} \
        -g {input.refined_bins_renamed}/*.fasta \
        -comp {params.completeness_min} \
        -sa {params.sa} \
        -nc {params.nc} \
        --genomeInfo {input.quality_summary_dRep_awk} \
        --checkM_method lineage_wf
        """

rule select_quality_bins:
# Among all the dereplicated MAGs, we'd like to know how many are medium to high-quality. See Bowers et al., 2017 Nature Biotechnology
    input:
        quality_summary = "results/{sample}/dereplicated_bins/data_tables/genomeInfo.csv"
    output:
        medium_quality_bins = "results/{sample}/medium_quality_bins.txt",
        high_quality_bins = "results/{sample}/high_quality_bins.txt",
        final_bin_set = "results/{sample}/final_bin_set.txt",
        final_bin_set_unique = "results/{sample}/final_bin_set_unique.txt"
    shell:
        """
        # Note that the 12th column is competeness and the 13th column is contamination. 
        awk -F "," '{{ if(($12 >= 50) && ($13 <10)) {{print}} }}' {input.quality_summary} | cut -f 1 > {output.medium_quality_bins}
        awk -F "," '{{ if(($12 > 90) && ($13 <5)) {{print}} }}' {input.quality_summary} | cut -f 1 > {output.high_quality_bins}
        cat {output.medium_quality_bins} {output.high_quality_bins} > {output.final_bin_set}
        sort {output.final_bin_set} | uniq > {output.final_bin_set_unique}
        """

rule copy_final_bins_over:
    input:
        bin_folder = "results/{sample}/dereplicated_bins/dereplicated_genomes/",
        final_bin_set_unique = "results/{sample}/final_bin_set_unique.txt"
    output:
        final_bin_folder = directory("results/{sample}/final_bin_set/"),
        copy_script = "results/{sample}/copy_final_bin_set.sh"
    shell:
        """
        mkdir -p {output.final_bin_folder}
        sed -e 's|^|cp {input.bin_folder}/|g' {input.final_bin_set_unique} > {output.copy_script}
        sed -i 's|$|.fasta {output.final_bin_folder}/.|g' {output.copy_script}
        bash {output.copy_script}
        """
        
rule assign_taxonomy:
# Now that we have the final bin set, let's assign their taxonomy.
    input:
        bins = "results/{sample}/final_bin_set/",
        final_bin_set_unique = "results/{sample}/final_bin_set_unique.txt"
    output:
        out_file = "results/{sample}/taxonomy.tsv",
        out_dir =  directory("results/{sample}/gtdbtk/")
    conda:
        "envs/gtdbtk-2.2.6.yml"
    params:
        ext = "fasta",
        threads = 20,
        db_path = "/storage1/data10/databases/release207_v2/"
    shell:
        """
        GTDBTK_DATA_PATH={params.db_path}

        gtdbtk classify_wf --skip_ani_screen \
        --genome_dir {input.bins} \
        --pplacer_cpus {params.threads} \
        -x {params.ext} --cpus {params.threads} \
        --out_dir {output.out_dir}

        grep -f {input.final_bin_set_unique} {output.out_dir}/gtdbtk/classify/*.summary.tsv >> {output.out_file}
        """