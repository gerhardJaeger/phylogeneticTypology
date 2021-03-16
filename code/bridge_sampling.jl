

function bridge_sampling(
      samples::Array{Float64,2},
      log_density::Function;
      verbose = false,
)
      nr = size(samples, 1) ÷ 2
      samples_4_fit = permutedims(Array(samples[1:nr, :]))
      samples_4_iter = permutedims(Array(samples[(nr+1):end, :]))
      r0 = 0.5
      tol1 = 1e-10
      n_post = size(samples_4_iter, 2)
      m = vec(mapslices(mean, samples_4_fit, dims = 2))::Vector{Float64}
      V = diagm(ones(size(samples, 2)))
      q11 = Vector{Float64}(undef, n_post) # unnormalized posterior
      @threads for i = 1:n_post
            q11[i] = log_density(Array(samples_4_iter[:, i]))
      end
      q12 = Vector{Float64}(undef, n_post) # g
      @threads for i = 1:n_post
            q12[i] = logpdf(MvNormal(m, V), samples_4_iter[:, i])
      end
      gen_samples = rand(MvNormal(m, V), n_post)
      q21 = Vector{Float64}(undef, n_post) # unnormalized posterior
      @threads for i = 1:n_post
            q21[i] = log_density(gen_samples[:, i])
      end
      q22 = Vector{Float64}(undef, n_post) # g
      @threads for i = 1:n_post
            q22[i] = logpdf(MvNormal(m, V), gen_samples[:, i])
      end
      l1 = q11 - q12 # samples_4_iter, unnormalized posterior / g
      l2 = q21 - q22# gen_samples, unnormalized posterior / g
      lstar = median(l1)
      n_1 = length(l1)
      n_2 = length(l2)
      neff = (1.0 * n_1)::Float64
      try
            neff = minimum(summarystats(sim).value[:, end, 1])
      catch e
      end
      s1 = neff ./ (neff .+ n_2)
      s2 = n_2 ./ (neff .+ n_2)
      logr = log(r0)
      log_r_vals = [logr]::Vector{Float64}
      logml = logr + lstar

      function upd(logr)
            lognumi = logsumexp([
                  l2[i] - lstar - logaddexp(log(s1) + l2[i] - lstar, log(s2) + logr) for
                  i = 1:size(l2, 1)
            ])
            logdeni = logsumexp([
                  -logaddexp(log(s1) + l1[i] - lstar, log(s2) .+ logr) for
                  i = 1:size(l1, 1)
            ])
            log(n_1) - log(n_2) + lognumi - logdeni
      end

      ft = optimize(y -> (upd(y)-y)^2, -10000., 0.)
      ft.minimizer + lstar

end
