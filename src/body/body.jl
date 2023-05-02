module body_module
"Define body strcuture and methods."

mutable struct body
  q::Vector{Float64}
  orientation
  q_old::Vector{Float64}
  orientation_old
  r_vectors::Array{Float64,2}
  Nmarkers::Int64

  function body(q, orientation, r_vectors)
    new(q, orientation, copy(q), copy(orientation), r_vectors, size(r_vectors)[1])
  end
end

function get_r_vectors(b)
  "Return r_vectors for the current body orientation."
  
  # Prepare vector and rotation matrix
  R = get_rotation_matrix(b.orientation)
  r_vectors = Array{Float64,2}(undef, 2, b.Nmarkers)

  # Compute points coordinates
  for i = 1 : b.Nmarkers
    r_vectors[1:end, i] = b.q + R * b.r_vectors[1:end, i]
  end

  return r_vectors
end

function get_rotation_matrix(orientation::Number)
  R = [cos(orientation) -sin(orientation);
       sin(orientation)  cos(orientation)]
  return R
end
                           


end


