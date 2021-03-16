using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
using Plots
using MCPhylo
using CSV
using DataFrames
using LinearAlgebra


##

function getEqulibrium(rates)
      q = zeros(4, 4)
      indices = [
            (2, 1),
            (3, 1),
            (4, 1),
            (1, 2),
            (3, 2),
            (4, 2),
            (1, 3),
            (2, 3),
            (4, 3),
            (1, 4),
            (2, 4),
            (3, 4),
      ]
      q[[CartesianIndex(x...) for x in indices]] = rates
      dia = sum(q, dims = 2)
      q[diagind(4, 4)] = -dia
      eq = svd(q).U[:, end]
      eq /= sum(eq)
end

##

fpairs = CSV.read("../data/fpairs.txt", DataFrame, header = false)[:, 1]

nrc = findall(occursin.("NRc", fpairs))
##



colnames = ["r_$(i)_$(j)" for j in 1:4 for i in 1:4 if i !=j]


##

df_ = []
for charNum = 1:28
      charName = fpairs[charNum]
      sim = read(
            "modelFitting/output/universal_$(lpad(charNum, 2, "0")).jls",
            ModelChains,
      )
      rates = exp.(vcat([sim[:, :s_rates, :].value[:, :, k] for k in sim.chains]...))
      if charNum ∈ nrc
            samples = samples[:, [4, 6, 5, 1, 3, 2, 11, 10, 12, 8, 7, 9]]
      end

      eqs = mapslices(getEqulibrium, rates, dims=2)

      rateMedians = mapslices(median, rates, dims=1)[1,:]
      rateHPD = permutedims(mapslices(hpd, rates, dims=1))

      rDF = DataFrame(
            :featurePair => charName,
            :variable => colnames,
            :median => rateMedians,
            :HPD_l => rateHPD[:,1],
            :HPD_u => rateHPD[:,2],
      )


      eqMedians = mapslices(median, eqs, dims=1)[1,:]
      eqHPD = permutedims(mapslices(hpd, eqs, dims=1))

      eqDF = DataFrame(
            :featurePair => charName,
            :variable => ["e1", "e2", "e3", "e4"],
            :median => eqMedians,
            :HPD_l => eqHPD[:,1],
            :HPD_u => eqHPD[:,2],
      )
      push!(df_, vcat(rDF, eqDF))
end

df = vcat(df_...)

CSV.write("../data/posteriorStatistics.csv", df)
