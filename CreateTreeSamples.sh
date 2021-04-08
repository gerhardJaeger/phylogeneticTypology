#!/bin/bash
cd code
# Download and create data fies
julia createData.jl
# run MrBayes for eligible files
# CAUTION: TIME INTENSIVE
julia runMrBayes.jl
# run RevBayes for for eligible files
# CAUTION: TIME INTENSIVE
cd revbayes
for f in *Rev; do
    rb $f;
done
# get posteriors from phylogenetic runs
cd ..
Rscript createPosterior.r
