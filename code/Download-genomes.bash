#!/usr/bin/env bash

echo "hello, I am downloading Ecoli genomes from NCBI"

#ncbi-genome-download --taxids 135618 bacteria # This is Methyloccales.
#cd data
ncbi-genome-download --taxids 511145 bacteria #This is ecoli.
#cd ..
