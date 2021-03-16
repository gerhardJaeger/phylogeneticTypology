## Code accompanying the paper *Phylogenetic typology* by Gerhard Jäger and Johannes Wahle



We used the following software:

- julia v.1.5.3
- R v.3.6.1 with libraries `ape`, `geiger`, `loo` and `LRO.utilities`
- MrBayes v. 3.2.7a with Beagle support
- RevBayes v.1.1.0

To replicate the results, change to the root directory of this repo and execute the following commands in that order:

- `cd code`

- `julia createData.jl`

- `cd ../data/asjpNex/`

- `for f in *mb.nex; do mb $f; done`

- `cd ../../code/revbayes`

- `for f in *Rev; do rb $f; done`

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

  

