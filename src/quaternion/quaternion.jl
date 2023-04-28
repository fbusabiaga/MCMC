# module quaternion_module
"Define quaternion and its methods."

# Load modules
import LinearAlgebra as LA


struct quaternion
  "Unit quaternion (s,p)."
  s::Float64
  p::Vector{Float64}

  function quaternion(s,p)
    q_norm = sqrt(s^2 + LA.dot(p,p))
    new(s / q_norm, p / q_norm)
  end

  function quaternion(q)
    q_norm = LA.norm(q)
    new(q[1] / q_norm, q[2:end] / q_norm)
  end 
end


function from_rotation(phi)
  "Create a quaternion given an angle of rotation phi"

  # Compute quaternion components
  phi_norm = LA.norm(phi)
  if phi_norm > 0
    s = cos(phi_norm / 2)
    p = (sin(phi_norm / 2) / phi_norm) * phi_norm
  else
    s = 1
    p = 0
  end

  # Create quaternion and return
  q = quaternion(s, p)
  return q
end


function rotation_matrix(q)
  "Return the rotation matrix representing rotation by this quaternion."

  # Get quaternion elements
  s = q.s
  p = q.p

  # Cross product matrix for p, actually the negative.
  diag = s^2 - 0.5
  R = 2.0 * [p[1]^2+diag       p[1]*p[2]-s*p[3]  p[1]*p[3]+s*p[2];
             p[2]*p[1]+s*p[3]  p[2]^2+diag       p[2]*p[3]-s*p[1];
             p[3]*p[1]-s*p[2]  p[3]*p[2]+s*p[1]  p[3]^2+diag]
  
  return R
end
  

# end
