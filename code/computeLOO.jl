global charNum = parse(Int, ARGS[1])

##
cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()


##
using LinearAlgebra
using StatsBase
using CSV
using DataFrames
using Pipe
using ProgressBars
using ProgressMeter
using StatsPlots
using Base.Threads
using RCall
using StatsFuns
using MCPhylo
using Random
Random.seed!(9011452010850473240)

##


wals = CSV.read("../data/charMtx.csv", DataFrame)
d = CSV.read("../data/fpairMtx.csv", DataFrame)

famFreqs = combine(groupby(wals, :glot_fam), nrow)
sort!(famFreqs, :nrow, rev=true)

families = filter(x -> x.nrow>1, famFreqs).glot_fam
isolates = filter(x -> x.nrow==1, famFreqs).glot_fam

isoDict = @pipe wals |>
      filter(x -> x.glot_fam ∈ isolates, _) |>
      zip(_.glot_fam, _.longname) |>
      Dict

##
fm2trees = Dict()
@showprogress for fm in families
      fm2trees[fm] = MCPhylo.ParseNewick("../data/posteriorTrees/$fm.posterior.tree")
end

for fm in isolates
      fm2trees[fm] = repeat([Node(isoDict[fm])], 1000)
end

##

function renumberTree(tree::GeneralNode)
      tips = [nd for nd in post_order(tree) if nd.nchild==0]
      nonTips = [nd for nd in post_order(tree) if nd.nchild > 0]
      sort!(tips, lt= (x,y) -> x.name < y.name)
      for (i,nd) in enumerate(tips)
            nd.num = i
      end
      for (i,nd) in enumerate(nonTips)
            nd.num = i+length(tips)
      end
      tree
end

##

for fm in keys(fm2trees)
      fm2trees[fm] = renumberTree.(fm2trees[fm])
end


lineages = vcat(families, isolates)


nLineages = length(lineages)


ttrees = [[fm2trees[fm][i] for fm in lineages] for i in 1:1000]



taxa = Vector{String}[]
for t in ttrees[1]
      push!(taxa, sort([nd.name for nd in pre_order(t) if nd.nchild==0]))
end


nnodes = [length(post_order(t)) for t in ttrees[1]]


states = ["a", "b", "c", "d"]

nsites = 1

nbase = 4

rates = [1.0]


##


data = zeros(nbase, nsites, maximum(nnodes), length(ttrees[1]))
char = Dict(zip(d.taxon, d[:,charNum+1]))
for i in 1:length(taxa)
      for l in taxa[i]
            s = char[l]
            t_ind = find_by_name(ttrees[1][i],l).num
            if s in states
                  data[:,1,t_ind, i] = (states .== s)
            else
                  data[:,1,t_ind, i] .= 1
            end
      end
end



##


uSim = read(
      "modelFitting/output/universal_$(lpad(charNum, 2, "0")).jls",
      ModelChains,
)


samples = vcat([uSim.value[:,:,k] for k in 1:2]...)

nSamples = size(samples, 1)

yPred = zeros(nSamples, nLineages)


@threads for i in ProgressBar(1:nSamples)
      s_rates = samples[i,2:(end-1)]
      ind = samples[i,1]
      treeI = ceil(Int, 1000*invlogit(ind))
      fullSrates = exp.(s_rates)
      y_ = zeros(nLineages)
      for fmI = 1:nLineages
            y_[fmI] = logpdf(
                  PhyloDist(ttrees[treeI][fmI], fullSrates, [1.0], freeK),
                  data[:, :, :, fmI],
            )
      end
      yPred[i,:] = y_

end

looU = convert(
      Dict,
      R"
      library(loo)
      loo($yPred)
      ",
)

##



lSim = read(
      "modelFitting/output/lineage_$(lpad(charNum, 2, "0")).jls",
      ModelChains,
)

samples = vcat([lSim.value[:,:,k] for k in 1:2]...)
samples = samples[sample(1:size(samples,1),1000),:]


nSamples = size(samples, 1)

yPred = zeros(nSamples, nLineages)


@threads for i = ProgressBar(1:nSamples)
      s_rates_ = samples[i,2:(end-1)]
      s_rates = reshape(s_rates_, :, nLineages)
      ind = samples[i,1]
      treeI = ceil(Int, 1000 * invlogit(ind))
      fullSrates = exp.(s_rates)
      y_ = zeros(nLineages)
      for fmI = 1:nLineages
            y_[fmI] = logpdf(
                  PhyloDist(ttrees[treeI][fmI], fullSrates[:, fmI], [1.0], freeK),
                  data[:, :, :, fmI],
            )
      end
      yPred[i, :] = y_
end

looL = convert(
      Dict,
      R"
      library(loo)
      loo($yPred)
      ",
)
##

kys = ["looic", "se_looic", "elpd_loo", "se_elpd_loo"]

results = permutedims(
      hcat(
            vcat(
                  [
                        [[k, looU[k], "u"] for k in kys],
                        [[k, looL[k], "l"] for k in kys],
                  ]...,
            )...,
      ),
)

resultsDF = DataFrame(results, :auto)
rename!(resultsDF, [:parameter, :value, :model])
resultsDF[:,:char] .= charNum

try mkdir("looFamResults")
catch e
end

CSV.write("looFamResults/$(lpad(charNum, 2, "0")).csv", resultsDF)
