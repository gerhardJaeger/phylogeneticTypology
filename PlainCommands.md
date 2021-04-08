# Plain Commands

This file lists all commands in order. These are the same commands as in the CreateTreeSamples.sh file and the ModelFitting.sh file.

- `cd code`

- `julia createData.jl`

- `julia runMrBayes.jl`

- `cd revbayes`

- `for f in *Rev; do rb $f; done`

- `cd ..`

- `Rscript createPosterior.r`

- `julia priorPredictiveSampling.jl`

- `cd modelFitting`

- `for i in {1..28}; do julia universal.jl $i; done`

- `for i in {1..28}; do julia lineage.jl $i; done`

- `cd ..`

- `for i in {1..28}; do julia posteriorPredictiveSampling.jl $i; done`

- `for i in {1..28}; do julia computeBS.jl $i; done`

- `for i in {1..28}; do julia computeLOO.jl $i; done`

- `julia correlations.jl`

- `julia visualizeBS.jl`

- `julia visualizeCTMC.jl`

- `julia visualizeLOO.jl`

- `julia visualizePPP.jl`

- `julia savePosteriorStatistics.jl`
