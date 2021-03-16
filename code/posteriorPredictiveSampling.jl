global charNum = parse(Int, ARGS[1])
##
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using CSV
using DataFrames
using Pipe
using Statistics
using ProgressMeter
using LinearAlgebra
using StatsFuns
using MCPhylo

##

wals = CSV.read("../data/charMtx.csv", DataFrame)

famFreqs = combine(groupby(wals, :glot_fam), nrow)
sort!(famFreqs, :nrow, rev=true)

families = filter(x -> x.nrow>1, famFreqs).glot_fam
isolates = filter(x -> x.nrow==1, famFreqs).glot_fam

lineages = vcat(families, isolates)

isoDict = @pipe wals |>
      filter(x -> x.glot_fam ∈ isolates, _) |>
      zip(_.glot_fam, _.longname) |>
      Dict

##

d = CSV.read("../data/fpairMtx.csv", DataFrame)
d = innerjoin(d, wals[:,[:longname, :glot_fam]], on=:taxon => :longname)
##

function myvar(x)
    if length(x) == 1
        return 0.
    end
    var(x, corrected=false)
end

##

fm2trees = Dict()
@showprogress for fm in families
      fm2trees[fm] = MCPhylo.ParseNewick("../data/posteriorTrees/$fm.posterior.tree")
end

for fm in isolates
      fm2trees[fm] = repeat([Node(isoDict[fm])], 1000)
end

##

offDiag = [CartesianIndex(x,y) for y in 1:4 for x in 1:4 if x!=y];
taxa = d.taxon;

defined = findall(d[:,charNum+1] .!= "-")

##

nSims = 1000

sim = read("modelFitting/output/universal_$(lpad(charNum, 2, "0")).jls", ModelChains)

samples = vcat([sim.value[:,:,k] for k in sim.chains]...)


##
stats = []
@showprogress for g in 1:nSims
    i = rand(1:size(samples, 1))
    treeI = ceil(Int,invlogit(samples[i,1])*1000)

    simChar = Dict{String,Int}()

    fullSrates = exp.(samples[i,2:(end-1)])
    q = zeros(4, 4)
    q[offDiag] = fullSrates
    for i = 1:4
        q[i, i] -= sum(q[i, :])
    end


    D, U = eigen(q)
    Uinv = inv(U)

    eq = real.(Uinv[4, :])
    eq /= sum(eq)

    for fmI = 1:length(lineages)
        fmI
        tree = fm2trees[lineages[fmI]][treeI]
        preo = pre_order(tree)
        tStates = Vector{Int}(undef, length(preo))
        tStates[tree.num] = rand(Categorical(eq))
        for nd in preo[2:end]
            mState = tStates[nd.mother.num]
            p = real.(U * diagm(exp.(nd.inc_length * D)) * Uinv)
            dState = rand(Categorical(p[mState, :]))
            tStates[nd.num] = dState
        end
        for nd in preo
            if nd.nchild == 0
                simChar[nd.name] = tStates[nd.num]
            end
        end
    end


    x = [simChar[l] for l in taxa]
    x_bin = DataFrame(hcat([x .== i for i = 1:4]...))
    x_bin[!, :glot_fam] = d.glot_fam
    x_bin = x_bin[defined,:]

    tv = mapslices(myvar, Array(x_bin[:, 1:4]), dims = 1) |> sum
    lmv = @pipe x_bin |>
          groupby(_, :glot_fam) |>
          combine(x -> sum(mapslices(myvar, Array(x[:, 1:4]), dims = 1)), _) |>
          mean(_[:, 2])

    vc = @pipe x_bin |>
          groupby(_, :glot_fam) |>
          combine(x -> mapslices(mean, Array(x[:, 1:4]), dims = 1), _) |>
          mapslices(myvar, Array(_[:, 2:end]), dims = 1) |>
          sum

    push!(stats, [tv, lmv, vc])
end

##
uSim = DataFrame(Tuple.(stats))
rename!(uSim, [:totalVariance, :lineagewiseVariance, :crossLineageVariance])

##

uSim[:,:method] .= "simulated"
uSim[:,:model] .= "universal"


##
sim = read("modelFitting/output/lineage_$(lpad(charNum, 2, "0")).jls", ModelChains)

samples = vcat([sim.value[:,:,k] for k in sim.chains]...)

##

stats = []
@showprogress for g in 1:nSims
    i = rand(1:size(samples, 1))
    treeI = ceil(Int,invlogit(samples[i,1])*1000)
    simChar = Dict{String,Int}()
    fullSrates = reshape(exp.(samples[i,2:(end-1)]), 12, :)


    for fmI in 1:length(lineages)
        q = zeros(4,4)
        q[offDiag] = fullSrates[:,fmI]
        for i in 1:4
            q[i,i] -= sum(q[i,:])
        end
        D, U = eigen(q)
        Uinv = inv(U)

        eq = abs.(Uinv[4,:])
        eq /= sum(eq)
        tree = fm2trees[lineages[fmI]][treeI]
        preo = pre_order(tree)
        tStates = Vector{Int}(undef, length(preo))
        tStates[tree.num] = rand(Categorical(eq/sum(eq)))
        for nd in preo[2:end]
              mState = tStates[nd.mother.num]
              p = real.(U * diagm(exp.(nd.inc_length * D)) * Uinv)
              ps = abs.(p[mState,:]/sum(p[mState,:]))
              ps /= sum(ps)
              dState = rand(Categorical(ps))
              tStates[nd.num] = dState
        end
        for nd in preo
            if nd.nchild == 0
                simChar[nd.name] = tStates[nd.num]
            end
        end
    end


    x = [simChar[l] for l in taxa]
    x_bin = DataFrame(hcat([x .== i for i in 1:4]...))
    x_bin[!,:glot_fam] = d.glot_fam
    x_bin = x_bin[defined,:]


    tv = mapslices(myvar, Array(x_bin[:,1:4]), dims=1) |> sum
    lmv = @pipe x_bin |> groupby(_, :glot_fam) |>
        combine(x -> sum(mapslices(myvar, Array(x[:,1:4]), dims=1)), _) |>
        mean(_[:,2])

    vc = @pipe x_bin |> groupby(_, :glot_fam) |>
        combine(x -> mapslices(mean, Array(x[:,1:4]), dims=1), _) |>
    mapslices(myvar, Array(_[:,2:end]), dims=1) |> sum

    push!(stats, [tv, lmv, vc])

end
##
lSim = DataFrame(Tuple.(stats))
rename!(lSim, names(uSim)[1:3])

lSim[:,:method] .= "simulated"
lSim[:,:model] .= "lineageSpecific"



##

try
    mkdir("ppp")
catch e
end

CSV.write("ppp/$(lpad(charNum,2,"0")).csv", vcat(uSim, lSim))
