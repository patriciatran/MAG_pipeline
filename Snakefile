rule download_ncbi_genomes:
    input:
        script = "code/Download-genomes.bash"
    output:
        directory("refseq") # Note to myself: The ncbi-download-genome script will create a folder in the current directory so we can either spell it out here or cd into the right one, or move it after the script runs.
    resources:
        tmpdir="/storage1/data10/tmp"
    conda:
        "envs/ncbi-genome-download.yml"
    shell:
        """
        TMPDIR="{resources.tmpdir}"
        {input.script}
        mv refseq/**/* data/.
        """