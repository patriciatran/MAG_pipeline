# About:
I built this snakemake pipeline to showcase how it can be used to take FASTQ files all the way into a set of dereplicated MAGs.

The general steps  are:
- Assembled the `fastq` reads using SPADES
- Bin the MAGs using `metawrap`
- Deplicate the MAGs (if relevant) using `dRep`.
- Determine MAG quality using `checkM`
- Select only MIMAG quality-standard MAGs for further analyses  (e.g. >50% complete, <10% contamination).
- Assign taxonomy of this MAG set using `gtdbtk`.

![DAG of the workflow](https://github.com/patriciatran/test_MAG_pipeline/blob/main/dag.svg)

# Status:
Work in progress

