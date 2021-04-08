using Pkg
Pkg.activate(".")
Pkg.instantiate()
##

using Plots
using MCPhylo
using CSV
using DataFrames
using LinearAlgebra

pyplot()
##

function getEqulibirium(rates)
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

theme(:solarized_light)
#
upscale = 2
fntsm = Plots.font("sans-serif", pointsize = round(10.0 * upscale))
fntlg = Plots.font("sans-serif", pointsize = round(15.0 * upscale))
default(titlefont = fntlg, guidefont = fntlg, tickfont = fntsm, legendfont = fntsm)
default(size = (1000 * upscale, 1000 * upscale)) #Plot canvas size
#sl_palette = PlotThemes._themes[:solarized_light].defaults[:palette]


##



sh = 0.05 # shift of arrows
rml = 1 # rate multiplier
brd = 0.3 # border around circles
dm = 30 * upscale # diameter

arPars = [2,1]


x = [0, 0, 1, 1]
y = [0, 1, 0, 1]

##
ctmc_ = []

for charNum = 1:28
      charName = fpairs[charNum]
      sim = read(
            "modelFitting/output/universal_$(lpad(charNum, 2, "0")).jls",
            ModelChains,
      )
      samples = vcat([sim[:, :s_rates, :].value[:, :, k] for k in sim.chains]...)
      if charNum ∈ nrc
            samples = samples[:, [4, 6, 5, 1, 3, 2, 11, 10, 12, 8, 7, 9]]
      end
      rates = exp.(mapslices(median, samples, dims = 1))[1, :]
      eqs = getEqulibirium(rates)
      features = string.(split(charName, "-"))
      xNames = [features[1][1:1], features[1][2:end]]
      yNames = [features[2][1:1], features[2][2:end]]
      p = plot(
            title = charName,
            legend = false,
            xlim = [-.2, 1.2],
            ylim = [-.2, 1.2],
            aspectratio = :equal,
            xticks = ([0, 1], [xNames[2] * xNames[1], xNames[1] * xNames[2]]),
            yticks = ([0, 1], [yNames[2] * yNames[1], yNames[1] * yNames[2]]),
            size=(4*500, 500*7),
      )
      for i = 1:4
            scatter!(x[i:i], y[i:i], markersize = dm * sqrt.(eqs[i:i]))
      end
      s, d = 2, 1
      plot!(
            [x[s], x[d]] .- sh,
            [y[s] - brd, y[d] + brd],
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[1],
      )
      s, d = 3, 1
      plot!(
            [x[s] - brd, x[d] + brd],
            [y[s], y[d]] .- sh,
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[2],
      )
      s, d = 4, 1
      plot!(
            [x[s] - brd / sqrt(2), x[d] + brd / sqrt(2)] .- sh / sqrt(2),
            [y[s] - brd / sqrt(2), y[d] + brd / sqrt(2)] .- -sh / sqrt(2),
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[3],
      )
      s, d = 1, 2
      plot!(
            [x[s], x[d]] .+ sh,
            [y[s] + brd, y[d] - brd],
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[4],
      )
      s, d = 3, 2
      plot!(
            [x[s] - brd / sqrt(2), x[d] + brd / sqrt(2)] .- sh / sqrt(2),
            [y[s] + brd / sqrt(2), y[d] - brd / sqrt(2)] .+ -sh / sqrt(2),
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[5],
      )
      s, d = 4, 2
      plot!(
            [x[s] - brd, x[d] + brd],
            [y[s], y[d]] .- sh,
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[6],
      )
      s, d = 1, 3
      plot!(
            [x[s] + brd, x[d] - brd],
            [y[s], y[d]] .+ sh,
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[7],
      )
      s, d = 2, 3
      plot!(
            [x[s] + brd / sqrt(2), x[d] - brd / sqrt(2)] .+ sh / sqrt(2),
            [y[s] - brd / sqrt(2), y[d] + brd / sqrt(2)] .- -sh / sqrt(2),
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[8],
      )
      s, d = 4, 3
      plot!(
            [x[s], x[d]] .- sh,
            [y[s] - brd, y[d] + brd],
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[9],
      )
      s, d = 1, 4
      plot!(
            [x[s] + brd / sqrt(2), x[d] - brd / sqrt(2)] .+ sh / sqrt(2),
            [y[s] + brd / sqrt(2), y[d] - brd / sqrt(2)] .+ -sh / sqrt(2),
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[10],
      )
      s, d = 2, 4
      plot!(
            [x[s] + brd, x[d] - brd],
            [y[s], y[d]] .+ sh,
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[11],
      )
      s, d = 3, 4
      plot!(
            [x[s], x[d]] .+ sh,
            [y[s] + brd, y[d] - brd],
            color = :black,
            arrow = arrow(arPars...),
            linewidth = rml * rates[12],
      )
      push!(ctmc_,p)
end
##
plot(ctmc_...,layout=grid(7,4, heights=repeat([1/7],28)))

savefig("../img/ctmc.pdf")
