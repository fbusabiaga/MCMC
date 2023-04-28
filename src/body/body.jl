module body_module
"Define body strcuture and methods."

# using TensorOperations

mutable struct body
  q::Vector{Float64}
  orientation
  r_vectors::Array{Float64,2}
end

function get_r_vectors(b)
  "Return r_vectors for the current body orientation."
  
  # Prepare vector and rotation matrix
  N_blobs = size(b.r_vectors, 1)
  R = get_rotation_matrix(b.orientation)
  r_vectors = Array{Float64,2}(undef, N_blobs, 2)
  
  # Compute points coordinates
  for i = 1 : N_blobs
    r_vectors[i,1:end] = b.q + R * b.r_vectors[i,1:end]
  end

  return r_vectors
end

function get_rotation_matrix(orientation::Float64)
  R = [cos(orientation) -sin(orientation);
       sin(orientation)  cos(orientation)]
  return R
end
                           


end


