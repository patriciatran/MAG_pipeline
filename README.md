# About:
I built this snakemake pipeline to showcase how FASTQ files can be taken all the way into a set of good-quality dereplicated MAGs.

![MAG pipeline logo](https://github.com/patriciatran/MAG_pipeline/blob/main/visuals/MAG_pipeline_logo.png)

The general steps  are:
- Assembled the `fastq` reads using `SPADES`
- Bin the MAGs using `metawrap`
- Refining MAGs using `dasTool`
- Deplicate the MAGs (if relevant) using `dRep`.
- Determine MAG quality using `checkM`
- Select only MIMAG quality-standard MAGs for further analyses  (e.g. >50% complete, <10% contamination).
- Assign taxonomy of this MAG set using `gtdbtk`.

![DAG of the workflow](https://github.com/patriciatran/test_MAG_pipeline/blob/main/dag.svg)

# Status:
Work in progress

# Thanks to:
This pipeline exists because of the folks making these programs available, please cite their work:
- SPADES: https://github.com/ablab/spades
- Metawrap: https://github.com/bxlab/metaWRAP
- Metabat1 and Metabat2: https://bitbucket.org/berkeleylab/metabat 
- Maxbin2: https://sourceforge.net/projects/maxbin/ 
- DasTool: https://github.com/cmks/DAS_Tool
- dRep: https://github.com/MrOlm/drep
- CheckM: https://github.com/Ecogenomics/CheckM 
- GTDBTK: https://github.com/Ecogenomics/GTDBTk

- Snakemake: https://snakemake.readthedocs.io/en/stable/
- conda/anaconda: https://docs.anaconda.com/anaconda/user-guide/faq/
- mamba: https://github.com/mamba-org/mamba

# Further Reading:
MIMAG Standards: https://www.nature.com/articles/nbt.3893
