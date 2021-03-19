#!/bin/bash
cd code
# Download and create data fies
julia createData.jl
# run MrBayes for eligible files
# CAUTION: TIME INTENSIVE
cd ../data/asjpNex/
for f in *mb.nex; do
    mb $f;
done
# run RevBayes for for eligible files
# CAUTION: TIME INTENSIVE
cd ../../code/revbayes
for f in *Rev; do
    rb $f;
done
# get posteriors from phylogenetic runs
cd ..
Rscript createPosterior.r
