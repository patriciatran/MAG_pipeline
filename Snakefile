# (c) Patricia Tran
# About: This is a snakemake pipeline that will take your FASTQ reads and turn them into a set medium and high quality dereplicated MAGs, with taxonomy.

### Workflow starts below ###

# Define input and output files
input_dir = "data/fastq"
output_dir = "results"

TMPDIR="/storage1/data10/tmp"

# from snakemake.io import glob_wildcards
# Double brackets in {{sample}} for subfolders
fastq_files = glob_wildcards(f"{input_dir}/{{sample}}_R1.fastq.gz")
samples = list(fastq_files.sample)

## Rules below ##
rule all:
    input:
        # from sympy import expand
        expand("results/{sample}/taxonomy.tsv", sample=samples)
        
# Define rules for each step in the workflow
rule spades:
# Description: This rules takes R1 and R2 paired FASTQ reads, assuming they are metagenomic, and assembles them using SPADES.s
    input:
        r1 = input_dir + "/{sample}_R1.fastq.gz",
        r2 = input_dir + "/{sample}_R2.fastq.gz",
    output:
        assembly = output_dir + "/{sample}/assembly"
    conda:
        "envs/spades.yml"
    params:
        memory = 100, #use max 600 on our server
        kmers = "33,55,77,99,127",
        threads = 15

    shell:
        """
        spades.py -m {params.memory} -k {params.kmers} \
        --meta -t {params.threads} -1 {input.r1} -2 {input.r2} \
        -o {output.assembly}
        """

rule unzip_reads:
    input: 
        r1 = input_dir + "/{sample}_R1.fastq.gz",
        r2 = input_dir + "/{sample}_R2.fastq.gz",
        assembly = output_dir + "/{sample}/assembly"
    output:
        r1_unzip =  input_dir + "/unzipped/{sample}_1.fastq",
        r2_unzip =  input_dir + "/unzipped/{sample}_2.fastq"
    shell:
        """
        gzip -c -d {input.r1} > {output.r1_unzip}
        gzip -c -d {input.r1} > {output.r2_unzip}
        """

rule metawrap:
# Using the assembled scaffods (reads) in the assembly,  and the reads, we will now bin them using metawrap, using 3 different softwares.
    input:
        assembly = "results/{sample}/assembly/contigs.fasta",
        r1_unzip =  input_dir + "/unzipped/{sample}_1.fastq",
        r2_unzip =  input_dir + "/unzipped/{sample}_2.fastq"
    output:
        bins = "results/{sample}/bins/"
    conda:
        "envs/metawrap_env_2.yml"
    params:
        threads = 15,
        tmpdir = "storage1/data10/tmp"

    shell:
        """
        TMPDIR={params.tmpdir}
        metawrap binning -o {output.bins} -a {input.assembly} \
        -t {params.threads} \
        --metabat1 --metabat2 --maxbin2 \
        {input.r1_unzip} {input.r2_unzip}
        """

rule refine_MAGS:
# Refine bins using dasTools. Requires to create a scaff to bin file first.
    input:
        bins = "results/{sample}/bins/",
        assembly = "results/{sample}/assembly/contigs.fasta"
    output:
        refined_bins_tables = "results/{sample}/refine_bins_tables/",
        refined_bins = "results/{sample}/refine_bins_DASTool_bins/"
    params:
        threads = 15
    conda:
        "envs/dasTool.yml"
    shell:
        """
        mkdir {output.refined_bins}
        Fasta_to_Scaffolds2Bin.sh -e fa -i {input.bins}/maxbin2_bins/ > {output.refined_bins_tables}/maxbin2_scaf2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -e fa -i {input.bins}/metabat1_bins/ > {output.refined_bins_tables}/metabat1_scaf2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -e fa -i {input.bins}/metabat2_bins/ > {output.refined_bins_tables}/metabat2_scaf2bin.tsv
        
        DAS_Tool --write_bins 1 --write_bin_evals 1 --create_plots 1 \
        -i {output.refined_bins_tables}/metabat1_scaf2bin.tsv,{output.refined_bins_tables}/metabat2_scaf2bin.tsv,{output.refined_bins_tables}/maxbin2_scaf2bin.tsv \
        -l metabat1,metabat2,maxbin2 \
        -c {input.assembly} \
        -o {output.refined_bins} \
        --threads {params.threads}
        """

rule rename_MAG:
# MAGs need to be called ".fasta" instead of ".fa" for dRep to work:
    input:
        refined_bins = "results/{sample}/refine_bins_DASTool_bins/"
    output:
        refined_bins_renamed = "results/{sample}/renamed_refined_MAGS/"
    shell:
        """
        for file in {input.refined_bins}/*.fa; do mv "$file" "${file%.fa}.fasta"; done        
        """

rule drep:
# Once we have the bins, we will dereplicate them as to not have repetitive genomes.
    input:
        refined_bins_renamed = "results/{sample}/renamed_refined_MAGS/"
    output:
        dereplicated_bins = "results/{sample}/dereplicated_bins/"
    conda:
        "envs/dRep.yml"
    params:
        completeness_min =  50, # Note that the defaut is 75
        sa = 0.99, # ANI threshold (default is 0.99)
        nc = 0.1, # coverage threshold (default is 0.1)
        tmp_dir = "/storage1/data10/tmp" #  Important because CheckM will eat your whole /tmp folder (default) otherwise
    shell:
        """
        TMPDIR={params.tmp_dir}
        dRep dereplicate \
        {output.dereplicated_bins} \
        -g {input.refined_bins_renamed}/*.fa \
        -comp {params.completeness_min} \
        -sa {params.sa} \
        -nc {params.nc} \
        -g {output.dereplicated_bins} \
        --checkM_method lineage_wf
        """

rule checkm:
# dRep technically already applies as CheckM quality check, but  we will run that on the dereplicated bins only.
    input:
        bins = "results/{sample}/dereplicated_bins/"
    output:
        folder = "results/{sample}/checkm/",
        quality_summary = "results/{sample}/checkm/quality_summary.tsv"
    conda:
        "envs/checkM.yml"
    params:
        threads = 15,
        extension = "fa", #file extension
        tmdir = "/storage1/data10/tmp",
        pplacer_threads = 15

    shell:
        """
        checkm lineage_wf -x {params.extension} \
         -t {params.threads} \
         -f {output.quality_summary} \
        {input.bins} {output.folder}
        """

rule select_quality_bins:
# Among all the dereplicated MAGs, we'd like to know how many are medium to high--quality. See Bowers et al., 2017 Nature Biotechnology
    input:
        quality_summary = "results/{sample}/checkm/quality_summary.tsv"
    output:
        medium_quality_bins = "results/{sample}/medium_quality_bins.txt",
        high_quality_bins = "results/{sample}/high_quality_bins.txt"
    shell:
        """
        # Note that the 12th column is competeness and the 13th column is contamination. 
        awk '$12 >= 50' {input.quality_summary} | awk '$13 < 10' | cut -f 1 > {output.medium_quality_bins}
        awk '$12 > 90' {input.quality_summary} | awk '$13 < 5' | cut -f 1 > {output.high_quality_bins}
        """
        

rule gtdbtk:
# Now that we have th
    input:
        medium_quality_bins = "results/{sample}/medium_quality_bins.txt",
        high_quality_bins = "results/{sample}/high_quality_bins.txt"
    output:
        taxonomy_medium_qual = "results/{sample}/taxonomy_medium_qual.tsv",
        taxonomy_high_qual = "results/{sample}/taxonomy_high_qual.tsv",
        out_file = "results/{sample}/taxonomy.tsv"
    conda:
        "envs/gtdbtk-1.7.0.yml"
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.medium_quality_bins} --out_dir {output.taxonomy_medium_qual}/gtdbtk --cpus 4 
        gtdbtk classify_wf --genome_dir {input.high_quality_bins} --out_dir {output.taxonomy_high_qual}/gtdbtk --cpus 4
        awk -F '\t' '$8 == \"High quality\" || $8 == \"Medium quality\"' {output.taxonomy_medium_qual}/gtdbtk/gtdbtk.ar122.classification_pplacer.tsv >> {output.out_file}
        """