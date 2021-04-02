cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
using MCPhylo
using CSV
using DataFrames
using Pipe

##
families = @pipe readdir("revbayes") |>
    split.(_, ".") |> first.(_) |> string.(_)
families = [x for x in families if x != "output"]
##
for fm in families
    c1 = select(
        CSV.read("revbayes/output/$(fm)_run_1.p", DataFrame),
        ["alpha", "pi[1]", "pi[2]", "rootAge"],
    )[501:end,:]
    c2 = select(
        CSV.read("revbayes/output/$(fm)_run_2.p", DataFrame),
        ["alpha", "pi[1]", "pi[2]", "rootAge"],
    )[501:end,:]
    sim = Chains(
        cat(Array(c1), Array(c2), dims=3),
        1,
    )
    gd = gelmandiag(sim)
    if maximum(gd.value[:,1]) > 1.1
        println("$fm did not converge")
    end
end
