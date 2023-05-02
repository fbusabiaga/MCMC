module many_bodies_utils
"Define utilities to handle many bodies."


function get_r_vectors(bodies)
  "Return r_vectors for the current bodies orientations."
 
  # Prepare vector and rotation matrix
  r_vectors = Main.body_module.get_r_vectors(bodies[1])
  offset = 1 + bodies[1].Nmarkers
  for b in bodies[2:end]
    r_vectors = hcat(r_vectors, Main.body_module.get_r_vectors(b))
    offset += b.Nmarkers
  end
  return r_vectors
end


function get_r_vectors(bodies, Nmarkers::Int64)
  "Return r_vectors for the current bodies orientations."

  # Prepare vector and rotation matrix
  r_vectors = Array{Float64,2}(undef, 2, Nmarkers)
  offset = 1
  for b in bodies
    r_vectors[1:end, offset:offset+b.Nmarkers-1] = body_module.get_r_vectors(b)
    offset += b.Nmarkers
  end
  
  return r_vectors
end


function get_r_vectors_offsets(bodies)
  "Return r_vectors offsets."
 
  # Prepare vector 
  r_vectors_offsets = Vector{Int64}([1])
  for b in bodies[2:end]
    push!(r_vectors_offsets, r_vectors_offsets[end] + b.Nmarkers)
  end
  return r_vectors_offsets
end


end

