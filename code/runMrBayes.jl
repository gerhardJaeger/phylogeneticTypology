cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

##
using CSV
using DataFrames
using Statistics
using Distributions
using Pipe


##
using Conda
Conda.pip_interop(true)
Conda.pip("install", "ete3")


ENV["PYTHON"] = ""
Pkg.build("PyCall")
using PyCall

ete3 = pyimport("ete3")


##
data = CSV.read("../data/charMtx.csv", DataFrame)

##
worldGlotF = download("https://osf.io/jyvgt/download", "../data/world_fullGlot.tre")

glot = ete3.Tree(worldGlotF)

glot.prune(data.longname)


##
families = open("../data/glot3.txt") do file
    readlines(file)
end


try
    mkdir("mrbayes/")
catch e
end

##

function mbScript(fm, ngen, append)
    fmTaxa = filter(x -> x.glot_fam == fm, data).longname
    nex = """
#Nexus
\tBegin MrBayes;
\t\tset seed=6789580436154794230;
\t\tset swapseed = 614090213;
\t\texecute ../data/asjpNex/$fm.nex;
\t\tlset rates=gamma coding=all;
"""

fmGlot = glot.copy()
fmGlot.prune(fmTaxa)
constraints = []
if length(fmTaxa) > 5
    for nd in fmGlot.get_descendants()
        if !nd.is_leaf()
            push!(constraints, nd.get_leaf_names())
        end
    end
end

    if length(constraints) > 0
        for (i, cn) in enumerate(constraints)
            nex *= "\t\tconstraint c$i = " * join(cn, " ") * ";\n"
        end

        nex *= "\t\tprset topologypr = constraints("
        nex *= join(["c$i" for i in 1:length(constraints)], ",") * ");\n"
    end

    nex *= """
\t\tprset brlenspr = clock:uniform;
\t\tprset clockvarpr = igr;
\t\tprset treeagepr=Gamma(0.05, 0.005);
\t\tprset shapepr=Exponential(10);
\t\tset beagleprecision=double;
\t\tmcmcp Burninfrac=0.5 stoprule=no stopval=0.01;
\t\tmcmcp filename=../data/asjpNex/output/$fm;
\t\tmcmcp samplefreq=1000 printfreq=5000 append=$append;
\t\tmcmc ngen=$ngen nchains=4 nruns=2;
\t\tsump;
\t\tsumt;
\tend;
"""
    nex
end


##
for fm in families
    mbFile = "mrbayes/$(fm).mb.nex"
    nrun = 1000000
    open(mbFile, "w") do file
        write(file, mbScript(fm, nrun, "no"))
    end
    command = `mpirun -np 8 mb $mbFile`
    run(command)
    function converged()
        pstat = CSV.read(
            "../data/asjpNex/output/$fm.pstat",
            DataFrame,
            header = 2,
            datarow = 3,
            delim = "\t",
            ignorerepeated = true,
        )

        maxPSRF = maximum(pstat.PSRF)

        tstat = CSV.read(
            "../data/asjpNex/output/$fm.tstat",
            DataFrame,
            header = 2,
            datarow = 3,
            delim = "\t",
            ignorerepeated = true,
        )

        meanStdev = mean(tstat[:,4])

        vstat = CSV.read(
            "../data/asjpNex/output/$fm.vstat",
            DataFrame,
            header = 2,
            datarow = 3,
            delim = "\t",
            ignorerepeated = true,
            missingstring="NA",
        ) |> dropmissing

        maxPSRF = maximum([maxPSRF, maximum(vstat.PSRF)])
        maxPSRF <= 1.1 && meanStdev <= 0.01
    end
    while !converged()
        nrun += 1000000
        open(mbFile, "w") do file
            write(file, mbScript(fm, nrun, "yes"))
        end
        run(command)
    end
end
