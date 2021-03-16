using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using DataFrames
using CSV
using Plots
using StatsPlots
using MCPhylo
using ProgressMeter
##

fpairs = CSV.read("../data/fpairs.txt", DataFrame, header=false)[:,1]

##

theme(:solarized_light)
#
upscale = 1 #8x upscaling in resolution
fntsm = Plots.font("sans-serif", pointsize=round(20.0*upscale))
fntlg = Plots.font("sans-serif", pointsize=round(25.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
default(size=(800*upscale,600*upscale)) #Plot canvas size

##

d = CSV.read("ppp/empiricalStats.csv", DataFrame)

d[:,:char] = 1:28
d[:,:method] .= "empirical"
d[:,:model] .= "nature"

for charNum in 1:28
    global d
    cdf = CSV.read("ppp/$(lpad(charNum, 2, "0")).csv", DataFrame)
    cdf[:,:char] .= charNum
    d = vcat(d, cdf)
end

##
n=28
##

tvPlots = []


for charNum in 1:n
    cDF = filter(x -> x.char == charNum, d)
    cl = filter(x -> x.model=="lineageSpecific", cDF).totalVariance
    cu = filter(x -> x.model=="universal", cDF).totalVariance
    eu = filter(x -> x.model=="nature", cDF).totalVariance
    p = plot([1,1], hpd(cu),
        size=(4*600, 400*7),
        linewidth=5,
        xlim=[0.5,2.5],
        legend=false,
        color=:black,
        xrotation=0,
        xticks=([1,2], ["", ""]),
        ylim=(0.27,0.85)
    )
    if charNum%4 != 2
        yticks!(p,[false])
    end
    plot!([1,1], hpd(cu, alpha=0.5), linewidth=15, color=:black, alpha=0.5)
    plot!([2,2], hpd(cl), linewidth=5, color=:black)
    plot!([2,2], hpd(cl, alpha=0.5), linewidth=15, color=:black, alpha=0.5)
    plot!(p,[.9,2.1],[eu[1],eu[1]], linewidth=2)
    title!(fpairs[charNum])
    push!(tvPlots, p)
end
##
xticks!(tvPlots[end-3], ([1,2], ["universal", "lineage-specific"]))
xticks!(tvPlots[end-2], ([1,2], ["universal", "lineage-specific"]))
xticks!(tvPlots[end-1], ([1,2], ["universal", "lineage-specific"]))
xticks!(tvPlots[end], ([1,2], ["universal", "lineage-specific"]))
p1 = plot(tvPlots...,
    layout=grid(7,4, heights=repeat([1/7],n)))

savefig("../data/img/pppTV.pdf")
##

lvPlots = []


for charNum in 1:n
    cDF = filter(x -> x.char == charNum, d)
    cl = filter(x -> x.model=="lineageSpecific", cDF).lineagewiseVariance
    cu = filter(x -> x.model=="universal", cDF).lineagewiseVariance
    eu = filter(x -> x.model=="nature", cDF).lineagewiseVariance
    p = plot([1,1], hpd(cu),
        size=(4*600, 400*7),
        linewidth=5,
        xlim=[0.5,2.5],
        legend=false,
        color=:black,
        xrotation=0,
        xticks=([1,2], ["", ""]),
        ylim=(0.03,0.16),
        )
    if charNum%4 != 2
        yticks!(p,[false])
    end
    plot!([1,1], hpd(cu, alpha=0.5), linewidth=15, color=:black, alpha=0.5)
    plot!([2,2], hpd(cl), linewidth=5, color=:black)
    plot!([2,2], hpd(cl, alpha=0.5), linewidth=15, color=:black, alpha=0.5)
    plot!(p,[.9,2.1],[eu[1],eu[1]], linewidth=2)
    title!(fpairs[charNum])
    push!(lvPlots, p)
end
##
xticks!(lvPlots[end-3], ([1,2], ["universal", "lineage-specific"]))
xticks!(lvPlots[end-2], ([1,2], ["universal", "lineage-specific"]))
xticks!(lvPlots[end-1], ([1,2], ["universal", "lineage-specific"]))
xticks!(lvPlots[end], ([1,2], ["universal", "lineage-specific"]))
p1 = plot(lvPlots...,
    layout=grid(7,4, heights=repeat([1/7],n)))


savefig("../data/img/pppLV.pdf")
##


clvPlots = []


for charNum = 1:n
    cDF = filter(x -> x.char == charNum, d)
    cl = filter(x -> x.model == "lineageSpecific", cDF).crossLineageVariance
    cu = filter(x -> x.model == "universal", cDF).crossLineageVariance
    eu = filter(x -> x.model == "nature", cDF).crossLineageVariance
    p = plot(
        [1, 1],
        hpd(cu),
        size=(4*600, 400*7),
        linewidth = 5,
        xlim = [0.5, 2.5],
        legend = false,
        color = :black,
        xrotation = 0,
        xticks = ([1, 2], ["", ""]),
        ylim = (0.28, 0.7),
    )
    if charNum % 4 != 2
        yticks!(p, [false])
    end
    plot!([1, 1], hpd(cu, alpha = 0.5), linewidth = 15, color = :black, alpha = 0.5)
    plot!([2, 2], hpd(cl), linewidth = 5, color = :black)
    plot!([2, 2], hpd(cl, alpha = 0.5), linewidth = 15, color = :black, alpha = 0.5)
    plot!(p, [0.9, 2.1], [eu[1], eu[1]], linewidth = 2)
    title!(fpairs[charNum])
    push!(clvPlots, p)
end
##
xticks!(clvPlots[end-3], ([1,2], ["universal", "lineage-specific"]))
xticks!(clvPlots[end-2], ([1,2], ["universal", "lineage-specific"]))
xticks!(clvPlots[end-1], ([1,2], ["universal", "lineage-specific"]))
xticks!(clvPlots[end], ([1,2], ["universal", "lineage-specific"]))
p1 = plot(clvPlots...,
    layout=grid(7,4, heights=repeat([1/7],n)))

savefig("../data/img/pppCLV.pdf")
##

function getStatistics(charNum, stats)
    df = select(filter(x -> x.char == charNum, d), [stats, :model])
    r = [filter(x -> x.model == "nature", df)[1,stats]]
    append!(r, hpd(filter(x -> x.model == "universal", df)[:,stats]))
    append!(r, hpd(filter(x -> x.model == "lineageSpecific", df)[:,stats]))
end

##

tv = DataFrame(Tuple.(getStatistics.(1:28, :totalVariance)))
rename!(
    tv,
    [:empirical, :universal_lower, :universal_upper, :lineage_lower, :lineage_upper],
)

insertcols!(tv, 1, :fpair => fpairs)

e = Array(rpad.(round.(tv[:,2], digits = 3), 5, "0"))

uHPD = mapslices(
    x -> "(" * join(x, ", ") * ")",
    Array(rpad.(round.(tv[:, 3:4], digits = 3), 5, "0")),
    dims = 2,
)

lHPD = mapslices(
    x -> "(" * join(x, ", ") * ")",
    Array(rpad.(round.(tv[:, 5:6], digits = 3), 5, "0")),
    dims = 2,
)

output = DataFrame(
    "feature pair" => tv.fpair,
    "empirical" => e,
    "universal model (HPD)" => uHPD[:,1],
    "lineage-specific model (HPD)" => lHPD[:,1],
)

CSV.write("../data/tables/tv.tex", output, newline="\\\\\n", delim="&")

##


lv = DataFrame(Tuple.(getStatistics.(1:28, :lineagewiseVariance)))
rename!(
    lv,
    [:empirical, :universal_lower, :universal_upper, :lineage_lower, :lineage_upper],
)

insertcols!(lv, 1, :fpair => fpairs)

e = Array(rpad.(round.(lv[:,2], digits = 3), 5, "0"))

uHPD = mapslices(
    x -> "(" * join(x, ", ") * ")",
    Array(rpad.(round.(lv[:, 3:4], digits = 3), 5, "0")),
    dims = 2,
)

lHPD = mapslices(
    x -> "(" * join(x, ", ") * ")",
    Array(rpad.(round.(lv[:, 5:6], digits = 3), 5, "0")),
    dims = 2,
)

output = DataFrame(
    "feature pair" => lv.fpair,
    "empirical" => e,
    "universal model (HPD)" => uHPD[:,1],
    "lineage-specific model (HPD)" => lHPD[:,1],
)

CSV.write("../data/tables/lv.tex", output, newline="\\\\\n", delim="&")
##


clv = DataFrame(Tuple.(getStatistics.(1:28, :crossLineageVariance)))
rename!(
    clv,
    [:empirical, :universal_lower, :universal_upper, :lineage_lower, :lineage_upper],
)

insertcols!(clv, 1, :fpair => fpairs)

e = Array(rpad.(round.(clv[:,2], digits = 3), 5, "0"))

uHPD = mapslices(
    x -> "(" * join(x, ", ") * ")",
    Array(rpad.(round.(clv[:, 3:4], digits = 3), 5, "0")),
    dims = 2,
)

lHPD = mapslices(
    x -> "(" * join(x, ", ") * ")",
    Array(rpad.(round.(clv[:, 5:6], digits = 3), 5, "0")),
    dims = 2,
)

output = DataFrame(
    "feature pair" => clv.fpair,
    "empirical" => e,
    "universal model (HPD)" => uHPD[:,1],
    "lineage-specific model (HPD)" => lHPD[:,1],
)

CSV.write("../data/tables/clv.tex", output, newline="\\\\\n", delim="&")

##
