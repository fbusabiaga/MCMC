"Define energy functions."

import LinearAlgebra as LA


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


function pairwise_energy(b, bodies, range)
  energy = 0
  for bi in bodies
    if bi === b
      continue
    end
    d = LA.norm(b.q - bi.q)
    denergy = d < range ? 1e+06 * (range - d) : 0
    energy += denergy
    # println("energy = ", d, " ", denergy)
  end
  return energy
end
