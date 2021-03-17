cd(@__DIR__)
using Pkg
Pkg.activate("..")
Pkg.instantiate()

##
using LinearAlgebra
using StatsBase
using CSV
using DataFrames
using Pipe
using ProgressMeter
using Random


##

using MCPhylo
Random.seed!(4928370335238343681)

##

wals = CSV.read("../../data/charMtx.csv", DataFrame)
d = CSV.read("../../data/fpairMtx.csv", DataFrame)

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
      fm2trees[fm] = MCPhylo.ParseNewick("../../data/posteriorTrees/$fm.posterior.tree")
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

##

lineages = vcat(families, isolates)

##

ttrees = [[fm2trees[fm][i] for fm in lineages] for i in 1:1000]



##

states = ["a", "b", "c", "d"]

nsites = 1

nbase = 4

rates = [1.0]

##
taxa = Vector{String}[]
for t in ttrees[1]
      push!(taxa, sort([nd.name for nd in pre_order(t) if nd.nchild==0]))
end

##
