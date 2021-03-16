global charNum = parse(Int, ARGS[1])
include("loadData.jl")

##

char = @pipe d |>
      zip(_[:, :taxon], _[:, charNum+1]) |>
      Dict
##



nnodes = [length(post_order(t)) for t in ttrees[1]]


data = zeros(nbase, nsites, maximum(nnodes), length(ttrees[1]))

for i in 1:length(taxa)
      for l in taxa[i]
            s = char[l]
            t_ind = find_by_name(ttrees[1][i],l).num
            if s in states
                  data[:,1,t_ind, i] = (states .== s)
            else
                  data[:,1,t_ind, i] .= 1
            end
      end
end

##


nLineages = length(lineages)


my_data = Dict{Symbol, Any}(
      :data => data,
      :nLineages => nLineages
)

##


model = Model(
      ind = Stochastic(0, () -> Logistic(), true),
      data = Stochastic(
            4,
            (fullSrates, nLineages, ind) -> MultiplePhyloDist(
                  ttrees[Integer(ceil(1000 * invlogit(ind[])))],
                  fullSrates,
                  ones(1, nLineages),
                  freeK,
            ),
            false,
      ),
      s_rates = Stochastic(
            2,
            nLineages -> UnivariateDistribution[
                  Normal(0,1) for i in 1:12, j in 1:nLineages
            ],
            true,
      ),
      fullSrates = Logical(
            2,
            s_rates -> exp.(s_rates),
            false,
      ),
)

##


inits = [Dict{Symbol, Any}(
      :data => data,
      :nLineages => nLineages,
      :s_rates => rand(Normal(0,1), (12, nLineages)),
      :ind => (a = zeros(); a[]=rand(Logistic()); a),
) for i in 1:2]

##



scheme = [Slice(:s_rates, 1.), Empirical(:ind, 1) ]

setsamplers!(model, scheme)

##

try
      mkdir("output")
catch e
end

##

sim = mcmc(
      model,
      my_data,
      inits,
      1000,
      burnin = 0,
      thin = 10,
      chains = 2,
      trees = false,
)
bi = 1 + size(sim)[1] ÷ 2
gd = gelmandiag(sim[bi:end,:,:])
psrf = maximum(gd.value[:,1])
write("output/lineage_$(lpad(charNum, 2, "0")).jls", sim[bi:end,:,:])
write("output/lineage_$(lpad(charNum, 2, "0")).log", string(psrf)*"\n")

while psrf > 1.1
      global sim, psrf, bi, gd
      @show psrf
      sim = mcmc(sim, 1000)
      bi = 1 + size(sim)[1] ÷ 2
      gd = gelmandiag(sim[bi:end,:,:])
      psrf = maximum(gd.value[:,1])
      write("output/lineage_$(lpad(charNum, 2, "0")).jls", sim[bi:end,:,:])
      write("output/lineage_$(lpad(charNum, 2, "0")).log", string(psrf)*"\n")
end
