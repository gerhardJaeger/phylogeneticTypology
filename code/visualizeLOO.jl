using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using DataFrames
using CSV
using Plots
using StatsPlots
using Plots.PlotMeasures
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

loo_ = []
for i in 1:28
    push!(loo_, CSV.read("looFamResults/$(lpad(i,2,"0")).csv", DataFrame))
end

looLong = vcat(loo_...)

loo = unstack(filter(x -> x.parameter == "elpd_loo", looLong), :char, :model, :value)

insertcols!(loo,1,:fpair=>fpairs)
select!(loo, [:fpair, :u, :l])
rename!(loo, [:fpair, :universal, :lineage])

loo[:,:delta] = loo.universal .- loo.lineage

sort!(loo, :delta)

##

scatter(loo.delta, 1:28,
    legend = false,
    yticks= (1:28,loo.fpair),
    xlab="Δ elpd",
    left_margin = upscale*5mm,
    markersize=upscale*3,
    )

savefig("../img/loo.pdf")

##

sort!(loo, :delta, rev=true)

output = DataFrame(
    "feature pair" => loo.fpair,
    "\$\\Delta\$ elpd" => round.(loo.delta, digits=1)
)

CSV.write("../data/tables/loo.tex", output, delim="&", newline="\\\\\n")
