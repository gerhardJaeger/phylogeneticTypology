#!/bin/bash
cd code

# Prior Predictive Sampling
julia priorPredictiveSampling.jl

cd modelFitting

# Fit Universal Model
for i in {1..28}; do
    julia universal.jl $i;
done

# Fit Lineage specific Model
# CAUTION: TIME INTENSIVE
for i in {1..28}; do
    julia lineage.jl $i;
done

# Posterior Predictive Sampling
cd ..

for i in {1..28}; do
    julia posteriorPredictiveSampling.jl $i;
done

# Bridge Sampling
for i in {1..28}; do
    julia computeBS.jl $i;
done

# LOO
for i in {1..28}; do
    julia computeLOO.jl $i;
done

# correlations
julia correlations.jl

# Visualization
julia visualizeBS.jl

julia visualizeCTMC.jl

julia visualizeLOO.jl

julia visualizePPP.jl

# save Results
julia savePosteriorStatistics.jl
