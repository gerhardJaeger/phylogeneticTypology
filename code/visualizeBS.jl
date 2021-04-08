using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using DataFrames
using CSV
using Plots
using StatsPlots
using Plots.PlotMeasures
using Statistics
using StatsFuns
##

fpairs = CSV.read("../data/fpairs.txt", DataFrame, header=false)[:,1]

##

theme(:solarized_light)
#
upscale = 2 #8x upscaling in resolution
fntsm = Plots.font("sans-serif", pointsize=round(8.0*upscale))
fntlg = Plots.font("sans-serif", pointsize=round(10.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
default(size=(400*upscale,600*upscale)) #Plot canvas size

##

bs_ = []
for i in 1:28
    push!(bs_, CSV.read("bsResults/$(lpad(i,2,"0")).csv", DataFrame, header=false))
end

bs = vcat(bs_...)
insertcols!(bs,1,:fpair=>fpairs)
rename!(bs, [:fpair, :universal, :lineage])

bs[:,:bf] = bs.universal .- bs.lineage

sort!(bs, :bf)

##

scatter(bs.bf, 1:28,
    legend = false,
    yticks= (1:28,bs.fpair),
    xlab="(log) Bayes factor",
    left_margin = upscale*5mm,
    markersize=upscale*3,
    )

savefig("../img/bayesfactor.pdf")



##

sort!(bs, :bf, rev=true)

p0 = logistic.(bs.bf)

bs[:,:cum_pv] = cumsum(1 .- p0)

output = DataFrame(
    "feature pair" => bs.fpair,
    "(log) Bayes Factor" => round.(bs.bf, digits=1),
    "cumulative posterior probability" => rpad.(round.(bs.cum_pv, sigdigits=3), 5, "0"),
)

CSV.write("../data/tables/bf.tex", output, delim="&", newline="\\\\\n")

##
