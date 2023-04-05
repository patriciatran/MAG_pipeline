rule all:
    input:
        directory("data/genomes"),
        "data/genome_names.txt"


rule download_ncbi_genomes:
    input:
        script = "code/Download-genomes.bash"
    output:
        directory("data/genomes")
        #"data/" + directory("genomes") # Note to myself: The ncbi-download-genome script will create a folder in the current directory so we can either spell it out here or cd into the right one, or move it after the script runs.
    resources:
        tmpdir="/storage1/data10/tmp"
    conda:
        "envs/ncbi-genome-download.yml"
    shell:
        """
        TMPDIR="{resources.tmpdir}"
        {input.script}
        mkdir data/genomes/
        mv refseq/bacteria/**/*.gz data/genomes/.
        """

rule get_genome_names:
    input:
        script = "code/Get_genome_names.bash"
    output:
        "data/genome_names.txt"
    params:

    shell:
        """
        {input.script} {output}
        """

# rule unzip_all_genomes:
#     input:
#         sample = "data/{sample}.gz"
#     output:
#         output = "data/unzipped/{sample}"
#     shell:
#         """
#         tar -xvg {input.sample}
#         """

# rule copyab:
#     input: 
#         "a.txt"
#     output: 
#         "b.txt"
#     shell:
#         """
#         copy {input} {output}
#         """

# rule copyac:
#     input:
#         "a.txt"
#     output:
#         # "c.txt"
#     shell:
#         "copy {input} {output}"