"Define energy functions."

import LinearAlgebra as LA


function external_enery(b, radius)
  r_vectors = body_module.get_r_vectors(b)
  r_norm = sqrt.(sum(abs2, r_vectors ,dims=1))
  r_max = maximum(r_norm)
  if r_max > radius
    return Inf
  else
    return 0
  end
end


function pairwise_energy_bodies(b, bodies, range)
  energy = 0
  for bi in bodies
    if bi === b
      continue
    end
    d = LA.norm(b.q - bi.q)
    denergy = d < range ? 1e+06 * (range - d) : 0
    energy += denergy
  end
  return energy
end


function pairwise_energy_surface(ri, r_vectors, e0, lambda)
  energy = 0
  for i = 1:size(ri, 2)
    rii = ri[1:end, i]
    for j = 1:size(r_vectors, 2)
      rj = r_vectors[1:end, j]
      d = LA.norm(rii - rj)
      energy += e0 * exp(-d / lambda)
    end
  end
  return energy
end
