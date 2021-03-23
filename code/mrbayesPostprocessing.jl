cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

##
using CSV
using DataFrames
using Statistics
using Distributions
using Pipe

##
pth = "../data/asjpNex/output"

pstatFiles = filter(x -> occursin(".pstat", x), readdir(pth, join=true))

#filter!(x -> !occursin("Austronesian.", x), pstatFiles)
#filter!(x -> !occursin("Atlantic-Congo.", x), pstatFiles)

pstat_ = []
for f in pstatFiles
    fm = first(split(last(split(f, "/")), "."))
    df = CSV.read(
        f,
        DataFrame,
        header = 2,
        datarow = 3,
        delim = "\t",
        ignorerepeated = true,
    )
    insertcols!(df, 1, :family=>fm)
    push!(pstat_, df)
end

pstat = vcat(pstat_...)

maxPSRF = combine(x -> maximum(x.PSRF), groupby(pstat, :family))

println(@pipe filter(x -> x.x1 > 1.1, maxPSRF).family |> string.(_))


##


tstatFiles = filter(x -> occursin(".tstat", x), readdir(pth, join = true))

#filter!(x -> !occursin("Austronesian.", x), tstatFiles)
#filter!(x -> !occursin("Atlantic-Congo.", x), tstatFiles)


tstat_ = []
for f in tstatFiles
    fm = first(split(last(split(f, "/")), "."))
    df = CSV.read(
        f,
        DataFrame,
        header = 2,
        datarow = 3,
        delim = "\t",
        ignorerepeated = true,
    )
    insertcols!(df, 1, :family => fm)
    push!(tstat_, df)
end

tstat = vcat(tstat_...)

meanStdv = combine(groupby(tstat[:, [1, 2, 5]], :family), "Stddev(s)" => mean)

println(meanStdv[meanStdv[:,2] .> 0.01,1])

##


vstatFiles = filter(x -> occursin(".vstat", x), readdir(pth, join = true))

vstat_ = []
for f in vstatFiles
    fm = first(split(last(split(f, "/")), "."))
    df = CSV.read(
        f,
        DataFrame,
        header = 2,
        datarow = 3,
        delim = "\t",
        ignorerepeated = true,
    )
    insertcols!(df, 1, :family => fm)
    push!(vstat_, df)
end

vstat = vcat(vstat_...)

println(string.(vstat[vstat.PSRF .> 1.1,:family] |> unique))
