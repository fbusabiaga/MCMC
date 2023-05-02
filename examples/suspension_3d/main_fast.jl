"Main code."

# Include external modules
import Random
import NearestNeighbors

# Include local modules
include("../../src/read_input/read_input.jl")
import .read_input
include("../../src/body/body.jl")
import .body_module
include("../../src/quaternion/quaternion.jl")
using .quaternion_module
include("../../src/many_bodies/many_bodies_utils.jl")
import .many_bodies_utils
include("energy.jl")

# Read input
options = read_input.read_input_file("data.main")

# Copy input file
output_name = haskey(options, "output_name") ? options["output_name"][1] : ".config" 
cp("data.main", output_name * ".inputfile", force=true)

# Init random numbers
rng = Random.MersenneTwister()

# Get MCMC options
kT = haskey(options, "kT") ? parse(Float64, options["kT"][1]) : 0
dx = haskey(options, "dx") ? parse(Float64, options["dx"][1]) : 1
dtheta = haskey(options, "dtheta") ? parse(Float64, options["dtheta"][1]) : pi

# Set energy parameters
radius = haskey(options, "radius") ? parse(Float64, options["radius"][1]) : Inf
e0 = haskey(options, "e0") ? parse(Float64, options["e0"][1]) : 0
lambda = haskey(options, "lambda") ? parse(Float64, options["lambda"][1]) : 1
d_tree = 2 * dx + 25 * lambda

# Create bodies
number_bodies = haskey(options, "number_bodies") ? parse(Int64, options["number_bodies"][1]) : 1
radius_body = haskey(options, "radius_body") ? parse(Int64, options["radius_body"][1]) : 1
r_vectors = transpose([0.000000000000000e+00	0.000000000000000e+00	1.000000000000000e+00;
                       2.763932022500211e-01	8.506508083520399e-01	4.472135954999580e-01;
                       -7.236067977499789e-01	5.257311121191337e-01	4.472135954999580e-01;
                       -7.236067977499790e-01	-5.257311121191336e-01	4.472135954999580e-01;
                       2.763932022500208e-01	-8.506508083520400e-01	4.472135954999580e-01;
                       8.944271909999159e-01	-2.190714793056811e-16	4.472135954999580e-01;
                       7.236067977499789e-01	5.257311121191337e-01	-4.472135954999580e-01;
                       -2.763932022500209e-01	8.506508083520400e-01	-4.472135954999580e-01;
                       -8.944271909999159e-01	1.095357396528405e-16	-4.472135954999580e-01;
                       -2.763932022500212e-01	-8.506508083520399e-01	-4.472135954999580e-01;
                       7.236067977499789e-01	-5.257311121191337e-01	-4.472135954999580e-01;
                       0.000000000000000e+00	0.000000000000000e+00	-1.000000000000000e+00])
bodies = []
for i = 1:number_bodies
  e_current = 1
  while e_current > 0
    global q = (Random.rand(rng, Float64, 3) - [0.5, 0.5, 0.5]) * 2 * radius
    e_current = LA.norm(q) > radius - radius_body ? Inf : 0
    for bi in bodies
      d = LA.norm(q - bi.q)
      denergy = d < (2 * radius_body + lambda) ? 1e+06 * ((2 * radius_body + lambda) - d) : 0
      e_current += denergy
    end
    println("i = ", i, ", e_current = ", e_current)
  end  
  orientation = quaternion_module.quaternion([1 0 0 0])  
  Nmarkers = 32
  theta = collect(range(0,2*pi,length=Nmarkers+1)[1:end-1])
  b = body_module.body(q, orientation, r_vectors)
  push!(bodies, b)
end

# Build r_vectors
r_vectors = many_bodies_utils.get_r_vectors(bodies)
r_vectors_offsets = many_bodies_utils.get_r_vectors_offsets(bodies)
Nmarkers = r_vectors_offsets[end]

# Open output files
f_config = open(output_name * ".config", "w")

# Loop steps
n_steps = parse(Int, options["n_steps"][1])
n_save = parse(Int, options["n_save"][1])
initial_step = haskey(options, "initial_step") ? parse(Int64, options["initial_step"][1]) : 0
for step = initial_step:1:n_steps-1
    
  # Save config
  if (step % n_save == 0) && (step >= 0)
    println("step = ", step)
    write(f_config, string(length(bodies)), "\n")
    for b in bodies
      println(f_config, join(b.q, " "), " ", b.orientation)
    end
  end

  # Build tree
  kdtree = NearestNeighbors.KDTree(r_vectors; leafsize = 64)

  # Loop over bodies
  bodies_indices = Random.shuffle(rng, Vector{Int64}(1:length(bodies)))
  for i = 1 : length(bodies)
    # Select body at random
    bi = bodies_indices[i]
    local b = bodies[bi]

    # Select its r_vectors
    ri = r_vectors[1:end, r_vectors_offsets[bi] : r_vectors_offsets[bi]+b.Nmarkers-1]

    # Find neighbors
    indx = NearestNeighbors.inrange(kdtree, ri, d_tree)

    # Compute energy
    e_current = external_enery(b, radius)
    e_current += pairwise_energy_surface(ri, r_vectors, indx, e0, lambda)
    
    # Random step
    dq = Random.randn(rng, Float64, 3) * dx
    b.q += dq
    dorientation = quaternion_module.quaternion(randn(rng, Float64, 4) * dtheta)
    b.orientation = dorientation * b.orientation    
    ri = body_module.get_r_vectors(b)
    r_vectors[1:end, r_vectors_offsets[bi] : r_vectors_offsets[bi]+b.Nmarkers-1] = body_module.get_r_vectors(b)

    # Compute new energy
    e_new = external_enery(b, radius)
    e_new += pairwise_energy_surface(ri, r_vectors, indx, e0, lambda)
  
    # Do Metropolies step
    mcmc = exp(-(e_new - e_current) / kT)
    rn = rand(rng, Float64)
    if rn < mcmc
      b.q_old = copy(b.q)
      b.orientation_old = deepcopy(b.orientation)
      r_vectors[1:end, r_vectors_offsets[bi] : r_vectors_offsets[bi]+b.Nmarkers-1] = body_module.get_r_vectors(b)
    else
      b.q = copy(b.q_old)
      b.orientation = deepcopy(b.orientation_old)
      r_vectors[1:end, r_vectors_offsets[bi] : r_vectors_offsets[bi]+b.Nmarkers-1] = body_module.get_r_vectors(b)
    end
  end
end

# Save last config
if (n_steps % n_save == 0)
  write(f_config, string(length(bodies)), "\n")
  for b in bodies
    println(f_config, join(b.q, " "), " ", b.orientation)
  end
end


