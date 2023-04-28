"Define energy functions."

function external_enery(b, radius)
  r_vectors = body_module.get_r_vectors(b)
  r_norm = sqrt.(sum(abs2, r_vectors ,dims=2))
  r_max = maximum(r_norm)

  if r_max > radius
    return Inf
  else
    return 0
  end
end
