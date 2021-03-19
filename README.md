## Code accompanying the paper *Phylogenetic typology* by Gerhard Jäger and Johannes Wahle


We used the following software:

- julia v.1.5.3
- R v.3.6.1 with libraries `ape`, `geiger`, `loo` and `LRO.utilities`
- MrBayes v. 3.2.7a with Beagle support
- RevBayes v.1.1.0

To replicate the results, change to the root directory of this repo and execute the following scripts:

- CreateTreeSamples.sh
- ModelFitting.sh

The first script will compute the posterior distributions of phylogenetic trees. The second script fits the
lineage specific and the universal model. It also compute the statistics for the model fit and produces
graphical and textual output.

For a plain list of the commands see here: [Plain Commands](./PlainCommands.md)
