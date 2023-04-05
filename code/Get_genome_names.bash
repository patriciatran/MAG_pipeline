#!/usr/bin/env bash

# This creates the column names
echo "genome_name" > data/genome_names.txt

# This uses tar to look at the file names, only get the name without the extension, and saves it to a file
ls -1 data/genomes/* | sed 's/.gbff.gz//g' | sed 's/.*G/G/g' >> data/genome_names.txt