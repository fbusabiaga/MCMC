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
    new(q, orientation, copy(q), deepcopy(orientation), r_vectors, size(r_vectors)[2])
  end
end

function get_r_vectors(b)
  "Return r_vectors for the current body orientation."
  
  # Prepare vector and rotation matrix
  R = get_rotation_matrix(b.orientation)
  r_vectors = Array{Float64,2}(undef, 2, b.Nmarkers)

  # Compute points coordinates
  r_vectors = b.q .+ R * b.r_vectors

  return r_vectors
end

function get_rotation_matrix(orientation::Number)
  R = [cos(orientation) -sin(orientation);
       sin(orientation)  cos(orientation)]
  return R
end

function get_rotation_matrix(orientation)
  R = Main.quaternion_module.rotation_matrix(orientation)
  return R
end



end


