charNum = parse(Int, ARGS[1])
charS = lpad(charNum, 2, "0")

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
using ProgressMeter
using StatsPlots
using MCPhylo
using Base.Threads
using StatsFuns
using Optim
using Distributions
using Random
Random.seed!(6897415724436294598)

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


function unnorm_log_post(mc::MCPhylo.ModelChains, θ_vec::Array)
      mod = deepcopy(mc.model)
      relist!(mod, θ_vec)
      logpdf(mod)
end


##

for fm in keys(fm2trees)
      fm2trees[fm] = renumberTree.(fm2trees[fm])
end

##

lineages = vcat(families, isolates)
nLineages = length(lineages)

##

ttrees = [[fm2trees[fm][i] for fm in lineages] for i in 1:1000]

nnodes = [length(post_order(t)) for t in ttrees[1]]

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

include("bridge_sampling.jl")
##

char = Dict(zip(d.taxon, d[:,charNum+1]))

data = zeros(nbase, nsites, maximum(nnodes), length(ttrees[1]))

for i = 1:length(taxa)
      for l in taxa[i]
            s = char[l]
            t_ind = find_by_name(ttrees[1][i], l).num
            if s in states
                  data[:, 1, t_ind, i] = (states .== s)
            else
                  data[:, 1, t_ind, i] .= 1
            end
      end
end

##
sim = read(
      "modelFitting/output/universal_$(charS).jls",
      ModelChains,
)
samples =
      vcat([sim[:, [:s_rates, :ind], :].value[:, :, k] for k = 1:2]...) |> DataFrame
function log_density(sampleRow)
      unnorm_log_post(sim, sampleRow)
end
uD = bridge_sampling(Array(samples), log_density)



##

sz =

sim = read(
      "modelFitting/output/lineage_$(charS).jls",
      ModelChains,
)
samples =
      vcat(
            [
                  sim[1:(size(sim.value)[1]÷1000):end, [:s_rates, :ind], :].value[
                        :,
                        :,
                        k,
                  ] for k = 1:2
            ]...,
      ) |> DataFrame
function log_density(sampleRow)
      unnorm_log_post(sim, sampleRow)
end
lD = bridge_sampling(Array(samples), log_density, verbose=true)

##


try
      mkdir("bsResults")
catch e
end

file = open("bsResults/$(charS).csv", "w")
write(file, join(string.([uD, lD]), ",")*"\n")
close(file)
