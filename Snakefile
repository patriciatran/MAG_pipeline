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
    input:
        r1 = input_dir + "/{sample}_R1.fastq.gz",
        r2 = input_dir + "/{sample}_R2.fastq.gz",
    output:
        assembly = output_dir + "/{sample}/assembly"
    conda:
        "envs/spades.yml"
    params:
        minimum_size = 2000,
        kmers = "33,55,77,99,127",
        threads = 5

    shell:
        "spades.py -m {params.minimum_size} -k {params.kmers}  --meta -t {params.threads} -1 {input.r1} -2 {input.r2} -o {output.assembly}"

rule metawrap:
    input:
        assembly = "results/{sample}/assembly",
        r1 = input_dir + "/{sample}_R1.fastq.gz",
        r2 = input_dir + "/{sample}_R2.fastq.gz"
    output:
        bins = "results/{sample}/bins/"
    conda:
        "envs/metawrap_env_2.yml"
    params:
        threads = 5

    shell:
        """
        metawrap binning -o {output.bins} -a {input.assembly} \
        -t {params.threads} \
        --metabat2 --metabat1 --maxbin2 \
        -1 {input.R1} -2 {input.R2}
        """

rule drep:
    input:
        bins = "results/{sample}/bins/*.fasta"
    output:
        dereplicated_bins = "results/{sample}/dereplicated_bins/"
    conda:
        "envs/dRep.yml"
    params:
        completeness_min =  50 # Note that the defaut is 75
        sa = 0.99 # ANI threshold (default is 0.99)
        nc = 0.1 # coverage threshold (default is 0.1)
        tmp_dir = "/storage1/data10/tmp" #  Important because CheckM will eat your whole /tmp folder (default) otherwise
    shell:
        """
        TMPDIR={params.tmp_dir}
        dRep dereplicate {input.bins} \
        -comp {params.completeness_min} \
        -sa {params.sa} \
        -nc {params.nc} \
        -g {output.dereplicated_bins} \
        --checkM_method lineage_wf
        """

rule checkm:
    input:
        bins = "results/{sample}/dereplicated_bins/"
    output:
        quality_summary = "results/{sample}/quality_summary.tsv"
    conda:
        "envs/checkM.yml"
    shell:
        "checkm lineage_wf -x fa -t 4 {input.bins} {output.quality_summary}"

rule select_quality_bins:
    input:
        quality_summary = "results/{sample}/quality_summary.tsv"
    output:
        medium_quality_bins = "results/{sample}/medium_quality_bins.txt",
        high_quality_bins = "results/{sample}/high_quality_bins.txt"
    shell:
        """
        awk '$12 >= 50' {input.quality_summary} | awk '$13 >= 90' | cut -f 1 > {output.medium_quality_bins}
        awk '$12 >= 75' {input.quality_summary} | awk '$13 >= 90' | cut -f 1 > {output.high_quality_bins}
        """
        

rule gtdbtk:
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