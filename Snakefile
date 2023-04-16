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
rule assemble_reads:
# Description: This rules takes R1 and R2 paired FASTQ reads, assuming they are metagenomic, and assembles them using SPADES.s
    input:
        r1 = input_dir + "/{sample}_R1.fastq.gz",
        r2 = input_dir + "/{sample}_R2.fastq.gz",
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

rule bin:
# Using the assembled scaffods (reads) in the assembly,  and the reads, we will now bin them using metawrap, using 3 different softwares.
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
        mkdir {output.refined_bins_tables}
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
        mkdir {output.refined_bins_renamed}
        rename 's/.fa/.fasta/' {input.refined_bins}/*.fa
        mv {input.refined_bins}/*.fasta {output.refined_bins_renamed}/.      
        """

rule dereplicate:
# Once we have the bins, we will dereplicate them as to not have repetitive genomes.
    input:
        refined_bins_renamed = "results/{sample}/renamed_refined_MAGS/"
    output:
        dereplicated_bins = directory("results/{sample}/dereplicated_bins/"),
        dereplicated_bins_folder = directory("results/{sample}/dereplicated_bins/dereplicated_genomes/")
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
        -g {input.refined_bins_renamed}/*.fasta \
        -comp {params.completeness_min} \
        -sa {params.sa} \
        -nc {params.nc} \
        --checkM_method lineage_wf
        """

rule quality_check:
# dRep technically already applies a CheckM quality check, but  we will run that on the dereplicated bins only.
    input:
        bins = "results/{sample}/dereplicated_bins/dereplicated_genomes/"
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
        {input.bins} {output.folder}
        """

rule select_quality_bins:
# Among all the dereplicated MAGs, we'd like to know how many are medium to high--quality. See Bowers et al., 2017 Nature Biotechnology
    input:
        quality_summary = "results/{sample}/checkm/quality_summary.tsv",
    output:
        medium_quality_bins = "results/{sample}/medium_quality_bins.txt",
        high_quality_bins = "results/{sample}/high_quality_bins.txt",
        final_bin_set = "results/{sample}/final_bin_set.txt",
        final_bin_set_unique = "results/{sample}/final_bin_set_unique.txt"
    shell:
        """
        # Note that the 12th column is competeness and the 13th column is contamination. 
        awk -F "\t" '{{ if(($12 >= 50) && ($13 <10)) {{print}} }}' {input.quality_summary} | cut -f 1 > {output.medium_quality_bins}
        awk -F "\t" '{{ if(($12 > 90) && ($13 <5)) {{print}} }}' {input.quality_summary} | cut -f 1 > {output.high_quality_bins}
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
        mkdir {output.final_bin_folder}
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

        gtdbtk classify_wf --skip-ani-screen \
        --genome_dir {input.bins} \
        --pplacer_cpus {params.threads} \
        -x {params.ext} --cpus {params.threads} \
        --out_dir {output.out_dir}

        grep -f {input.final_bin_set_unique} {output.out_dir}/gtdbtk/classify/*.summary.tsv >> {output.out_file}
        """