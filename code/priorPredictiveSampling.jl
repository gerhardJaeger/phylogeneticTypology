using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using CSV
using DataFrames
using Plots
using StatsPlots
using Pipe
using Statistics
using RCall
using ProgressMeter
using LinearAlgebra
using StatsFuns
using MCPhylo
using Random
Random.seed!(8074252001070266643)

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



function getEmpiricalStatistics(charNum)
    fp = names(d)[charNum+1]
    f1, f2 = Symbol.(split(fp, "-"))
    binaryRepresentation = select(
        innerjoin(
            d[:, [1, charNum + 1]],
            wals[:, [:glot_fam, :longname, f1, f2]],
            on = :taxon => :longname,
        ) |> dropmissing,
        ["taxon", "glot_fam", fp],
    )
    binaryRepresentation[!, :a] = Int.(binaryRepresentation[:, 3] .== "a")
    binaryRepresentation[!, :b] = Int.(binaryRepresentation[:, 3] .== "b")
    binaryRepresentation[!, :c] = Int.(binaryRepresentation[:, 3] .== "c")
    binaryRepresentation[!, :d] = Int.(binaryRepresentation[:, 3] .== "d")
    totalVariance = sum(mapslices(myvar, Array(binaryRepresentation[:, 4:end]), dims = 1))
    lineagewiseVariance = combine(
        x -> sum(mapslices(myvar, Array(x[:, 4:end]), dims = 1)),
        groupby(binaryRepresentation, :glot_fam),
    )
    rename!(lineagewiseVariance, [:glot_fam, :variance])
    lvar = mean(lineagewiseVariance.variance)
    varianceBetweenLineages = @pipe binaryRepresentation |>
          groupby(_, :glot_fam) |>
          combine(x -> mapslices(mean, Array(x[:, 4:end]), dims = 1), _) |>
          Array(_[:, 2:end]) |>
          mapslices(myvar, _, dims = 1) |>
          sum
    totalVariance, lvar, varianceBetweenLineages
end

empiricalStats = DataFrame(getEmpiricalStatistics.(1:28))
rename!(empiricalStats, [:totalVariance, :lineagewiseVariance, :crossLineageVariance])


##

treefiles = readdir("../data/posteriorTrees", join=true);

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

##

nSims = 1000

stats = []

##
@showprogress for g in 1:nSims
    treeI = sample(1:1000)

    simChar = Dict{String,Int}()

    global fullSrates = rand(LogNormal(0, 1), 12)
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
rename!(uSim, names(empiricalStats))

##

df = copy(empiricalStats)
df[:,:method] .= "empirical"
df[:,:model] .= "nature"

uSim[:,:method] .= "simulated"
uSim[:,:model] .= "universal"

df = vcat(df, uSim)

##

nSims = 1000

lstats = zeros(nSims, 3)

@showprogress for g in 1:nSims
    global q, fullSrates, p
    treeI = sample(1:1000)
    simChar = Dict{String,Int}()
    for fmI in 1:length(lineages)
        fullSrates = rand(LogNormal(0,1), 12)
        q = zeros(4,4)
        q[offDiag] = fullSrates
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
    x_bin

    tv = mapslices(myvar, Array(x_bin[:,1:4]), dims=1) |> sum
    lmv = @pipe x_bin |> groupby(_, :glot_fam) |>
        combine(x -> sum(mapslices(myvar, Array(x[:,1:4]), dims=1)), _) |>
        mean(_[:,2])

    vc = @pipe x_bin |> groupby(_, :glot_fam) |>
        combine(x -> mapslices(mean, Array(x[:,1:4]), dims=1), _) |>
    mapslices(myvar, Array(_[:,2:end]), dims=1) |> sum

    lstats[g,:] = [tv, lmv, vc]
end
##
lSim = DataFrame(lstats)
rename!(lSim, names(empiricalStats))

lSim[:,:method] .= "simulated"
lSim[:,:model] .= "lineageSpecific"

df = vcat(df, lSim)

##

theme(:solarized_light)

upscale = 1 #8x upscaling in resolution
fntsm = Plots.font("sans-serif", pointsize=round(20.0*upscale))
fntlg = Plots.font("sans-serif", pointsize=round(25.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
default(size=(800*upscale,600*upscale)) #Plot canvas size

##

u_tv_df = @pipe df |>
      filter(x -> x.model != "lineageSpecific", _) |>
      select(_, [:totalVariance, :method]) |>
      stack

p1 = @df u_tv_df boxplot(
    :method,
    :value,
    legend = false,
    grouping = :method,
    title = "total variance",
    linewidth=2*upscale
)

@df u_tv_df dotplot!(p1, :method, :value, markersize=2*upscale, color=:black, alpha=0.4)


##

u_lv_df = @pipe df |>
      filter(x -> x.model != "lineageSpecific", _) |>
      select(_, [:lineagewiseVariance, :method]) |>
      stack

p2 = @df u_lv_df boxplot(
    :method,
    :value,
    legend = false,
    grouping = :method,
    title = "lineage-wise variance",
    linewidth=2*upscale,
)
@df u_lv_df dotplot!(p2, :method, :value, markersize=2*upscale, color=:black, alpha=0.4)

##

u_clv_df = @pipe df |>
      filter(x -> x.model != "lineageSpecific", _) |>
      select(_, [:crossLineageVariance, :method]) |>
      stack

p3 = @df u_clv_df boxplot(
    :method,
    :value,
    legend = false,
    grouping = :method,
    title = "cross-lineage variance",
    linewidth=2*upscale
)
@df u_clv_df dotplot!(p3, :method, :value, markersize=2*upscale, color=:black, alpha=0.4)

##

savefig(p1, "../../img/priorSimulationUniversal_1.png")
savefig(p2, "../../img/priorSimulationUniversal_2.png")
savefig(p3, "../../img/priorSimulationUniversal_3.png")

##


l_tv_df = @pipe df |>
      filter(x -> x.model != "universal", _) |>
      select(_, [:totalVariance, :method]) |>
      stack

p4 = @df l_tv_df boxplot(
    :method,
    :value,
    legend = false,
    grouping = :method,
    title = "total variance",
    linewidth=2*upscale
)

@df l_tv_df dotplot!(p4, :method, :value, markersize=2*upscale, color=:black, alpha=0.4)


##

l_lv_df = @pipe df |>
      filter(x -> x.model != "universal", _) |>
      select(_, [:lineagewiseVariance, :method]) |>
      stack

p5 = @df l_lv_df boxplot(
    :method,
    :value,
    legend = false,
    grouping = :method,
    title = "lineage-wise variance",
    linewidth=2*upscale,
)
@df l_lv_df dotplot!(p5, :method, :value, markersize=2*upscale, color=:black, alpha=0.4)

##

l_clv_df = @pipe df |>
      filter(x -> x.model != "universal", _) |>
      select(_, [:crossLineageVariance, :method]) |>
      stack

p6 = @df l_clv_df boxplot(
    :method,
    :value,
    legend = false,
    grouping = :method,
    title = "cross-lineage variance",
    linewidth=2*upscale
)
@df l_clv_df dotplot!(p6, :method, :value, markersize=2*upscale, color=:black, alpha=0.4)

##

savefig(p4, "../../img/priorSimulationLineage_1.png")
savefig(p5, "../../img/priorSimulationLineage_2.png")
savefig(p6, "../../img/priorSimulationLineage_3.png")
