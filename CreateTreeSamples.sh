#!/bin/bash
cd code
# Download and create data fies
julia createData.jl
# run MrBayes for eligible files
# CAUTION: TIME INTENSIVE
cd mrbayes
for f in *mb.nex; do
    bn=`basename $f .mb.nex`;
    mb $f;
    mb $bn".mb1.nex";
done
# run RevBayes for for eligible files
# CAUTION: TIME INTENSIVE
cd ../revbayes
for f in *Rev; do
    rb $f;
done
# get posteriors from phylogenetic runs
cd ..
Rscript createPosterior.r
