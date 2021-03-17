
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using DataFrames
using CSV
using Plots
using StatsPlots
using Plots.PlotMeasures
using MCPhylo
using LinearAlgebra
using RCall
using Formatting
using ProgressMeter
using Statistics
using Random
Random.seed!(818480685100856001)

##

fpairs = CSV.read("../data/fpairs.txt", DataFrame, header=false)[:,1]

##
function getEqulibirium(rates)
      q = zeros(4, 4)
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
      eq = svd(q).U[:,end]
      eq /= sum(eq)
end

##

function cor4(eq)
      cv = eq[1]*eq[4] - eq[2]*eq[3]
      vx = (eq[1]+eq[2])*(eq[3]+eq[4])
      vy = (eq[1]+eq[3])*(eq[2]+eq[4])
      cv/sqrt(vx*vy)
end

##
nrc = findall(occursin.("NRc",fpairs))

##
eqs_ = []

for charNum = 1:28
      sim = read("modelFitting/output/universal_$(lpad(charNum,2,"0")).jls", ModelChains)
      s_rates = vcat([sim[:, :s_rates, :].value[:, :, k] for k = 1:2]...)
      posteriorEq = mapslices(getEqulibirium, exp.(s_rates), dims = 2)
      if charNum ∈ nrc
            posteriorEq = posteriorEq[:, [2, 1, 4, 3]]
      end
      posteriorCor = mapslices(cor4, posteriorEq, dims = 2)[:, 1]
      push!(eqs_, mapslices(mean, posteriorEq, dims = 1)[1, :])
end

eqs = DataFrame(Tuple.(eqs_))
rename!(eqs, ["00", "01", "10", "11"])
insertcols!(eqs, 1, :fpair => fpairs)



##

cors_ = []

for charNum = 1:28
      sim = read("modelFitting/output/universal_$(lpad(charNum,2,"0")).jls", ModelChains)
      s_rates = vcat([sim[:, :s_rates, :].value[:, :, k] for k = 1:2]...)
      posteriorEq = mapslices(getEqulibirium, exp.(s_rates), dims = 2)
      if charNum ∈ nrc
            posteriorEq = posteriorEq[:, [2, 1, 4, 3]]
      end
      posteriorCor = mapslices(cor4, posteriorEq, dims = 2)[:, 1]
      push!(
            cors_,
            vcat(
                  hpd(posteriorCor, alpha=0.05),
                  hpd(posteriorCor, alpha = 0.5),
                  [median(posteriorCor)],
            ),
      )
end


cors = DataFrame(Tuple.(cors_))
rename!(cors, [:hpd95_l, :hpd95_u, :hpd50_l, :hpd50_u, :median])
insertcols!(cors, 1, :fpair => fpairs)

sort!(cors, :median)

##


theme(:solarized_light)
#
upscale = 2 #8x upscaling in resolution
fntsm = Plots.font("sans-serif", pointsize=round(15.0*upscale))
fntlg = Plots.font("sans-serif", pointsize=round(20.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
default(size=(800*upscale,1000*upscale)) #Plot canvas size
sl_palette = PlotThemes._themes[:solarized_light].defaults[:palette]

##

linreg_ = []

for (g,charNum) in enumerate(indexin(cors.fpair, fpairs))
      charName = fpairs[charNum]
      features = string.(split(charName, "-"))
      xNames = [features[1][1:1], features[1][2:end]]
      yNames = [features[2][1:1], features[2][2:end]]
      sim = read("modelFitting/output/universal_$(lpad(charNum,2,"0")).jls", ModelChains)
      s_rates = vcat([sim[:, :s_rates, :].value[:, :, k] for k = 1:2]...)
      posteriorEq = mapslices(getEqulibirium, exp.(s_rates), dims = 2)
      if charNum ∈ nrc
            posteriorEq = posteriorEq[:, [2, 1, 4, 3]]
      end
      ms = 150 * upscale
      p = plot(
            aspectratio=:equal,
            legend = false,
            xlim = [-0.5, 1.5],
            ylim = [-0.5, 1.5],
            xticks = ([0, 1], [xNames[2] * xNames[1], xNames[1] * xNames[2]]),
            yticks = ([0, 1], [yNames[2] * yNames[1], yNames[1] * yNames[2]]),
            left_margin = upscale * 40mm,
            right_margin = upscale * 20mm,
            size=(4*500, 500*7),
      )
#      title!(fpairs[charNum])
      for g in 1:100
            i = rand(1:size(posteriorEq, 1))
            eqs = posteriorEq[i, :]
            scatter!(
                  [0],
                  [0],
                  markersize = ms * sqrt(eqs[1]),
                  alpha = 0.05,
                  markerstrokewidth = 0,
                  color = sl_palette.colors[1],
            )
            scatter!(
                  [0],
                  [1],
                  markersize = ms * sqrt(eqs[2]),
                  alpha = 0.05,
                  markerstrokewidth = 0,
                  color = sl_palette.colors[2],
            )
            scatter!(
                  [1],
                  [0],
                  markersize = ms * sqrt(eqs[3]),
                  alpha = 0.05,
                  markerstrokewidth = 0,
                  color = sl_palette.colors[3],
            )
            scatter!(
                  [1],
                  [1],
                  markersize = ms * sqrt(eqs[4]),
                  alpha = 0.05,
                  markerstrokewidth = 0,
                  color = sl_palette.colors[4],
            )
            a = eqs[2]/(eqs[1]+eqs[2])
            b = eqs[4]/(eqs[3]+eqs[4]) - a
            Plots.abline!(b, a, color=:gray, linewidth=upscale=0.5)
            end
      push!(linreg_, p)
end
plot(linreg_...,layout=grid(7,4, heights=repeat([1/7],28)))

try
      mkdir("../data/img/")
catch e
end


savefig("../data/img/linearRegression.pdf")
##

sort!(cors, :median)
p = plot(legend = false, left_margin = upscale*10mm)
for charNum = 1:28
      plot!(
            p,
            Array(cors[charNum, 2:3]),
            [charNum, charNum],
            linewidth = 3 * upscale,
            color = :black,
      )
      plot!(
            p,
            Array(cors[charNum, 4:5]),
            [charNum, charNum],
            linewidth = 10 * upscale,
            color = :black,
      )
      scatter!(
            p,
            [cors[charNum, 6]],
            [charNum],
            markersize = 5 * upscale,
            color = :white,
      )
end

yticks!(p, (1:28, cors.fpair))
vline!(p,[0])
xlabel!("correlation (posterior distribution)")

display(p)

savefig("../data/img/correlations.pdf")

##

sort!(cors, :median, rev=true)

output = DataFrame(
      "feature pair" => cors.fpair,
      "median" => format.("{1:.2f}",cors.median),
      "HPD" => mapslices(
            x -> "(" * join(x, ", ") * ")",
            format.("{1:.2f}", Array(cors[:, 2:3])),
            dims = 2,
      )[:,1],
)

try
      mkdir("../data/tables/")
catch e
end

CSV.write("../data/tables/correlations.tex", output, newline="\\\\\n", delim="&")


##

prior = rand(Normal(0,1), (100000, 12))

priorEqs = mapslices(getEqulibirium, exp.(prior), dims=2)

priorCor = mapslices(cor4, priorEqs, dims=2)[:,1]

##
bf_ = []

for charNum = 1:28
      sim = read(
            "modelFitting/output/universal_$(lpad(charNum, 2, "0")).jls",
            ModelChains,
      )
      s_rates = vcat([sim[:, :s_rates, :].value[:, :, k] for k = 1:2]...)
      posteriorEq = mapslices(getEqulibirium, exp.(s_rates), dims = 2)
      if charNum ∈ nrc
            posteriorEq = posteriorEq[:, [2, 1, 4, 3]]
      end
      posteriorCor = mapslices(cor4, posteriorEq, dims = 2)[:, 1]
      bfc = convert(
            Float64,
            R"library('LRO.utilities')
            log(savage_dickey($posteriorCor, $priorCor, Q=0, plot=F)$BF01)",
      )
      push!(bf_, bfc)
end
##

bf = DataFrame(:fpair => fpairs, :bf => bf_)

sort!(bf, :bf, rev = true)


p0 = invlogit.(bf.bf)

bf[:, :cum_pv] = [minimum([1.,x]) for x in cumsum(1 .- p0)]

output = DataFrame(
      "feature pair" => bf.fpair,
      "(log) Bayes Factor" => format.("{1:.2f}",round.(bf.bf, digits = 2)),
      "cumulative posterior probability" =>
            rpad.(round.(bf.cum_pv, sigdigits = 3), 5, "0"),
)

CSV.write("../data/tables/bfCorr.tex", output, delim="&", newline="\\\\\n")

##
upscale = 2

sort!(bf, :bf)

scatter(bf.bf, 1:28,
    legend = false,
    yticks= (1:28,bf.fpair),
    xlab="(log) Bayes factor",
    left_margin = upscale*5mm,
    markersize=upscale*3,
    )

savefig("../data/img/bayesfactorCorr.pdf")
