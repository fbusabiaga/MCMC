"Main code."

# Include external modules
import Random

# Include local modules
include("../../src/read_input/read_input.jl")
import .read_input
include("../../src/body/body.jl")
using .body_module
include("../../src/quaternion/quaternion.jl")
using .quaternion_module
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
radius = haskey(options, "radius") ? parse(Int, options["radius"][1]) : Inf

# Create body
q = [0, 0, 0]
orientation = quaternion_module.quaternion([1 0 0 0])
N_blobs = 32
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
b = body_module.body(q, orientation, r_vectors)

# Open output files
f_config = open(output_name * ".config", "w")

# Loop steps
n_steps = parse(Int, options["n_steps"][1])
n_save = parse(Int, options["n_save"][1])
initial_step = haskey(options, "initial_step") ? parse(Int64, options["initial_step"][1]) : 0
for step = initial_step:1:n_steps
  println("step = ", step)

  # Save config
  if (step % n_save == 0) && (step >= 0)
    write(f_config, string(1), "\n")
    println(f_config, join(b.q, " "), " ", b.orientation)
  end
  
  # Compute energy
  e_current = external_enery(b, radius)

  # Random step
  dq = randn(rng, Float64, 3) * dx
  b.q += dq
  dorientation = quaternion_module.quaternion(randn(rng, Float64, 4) * dtheta)
  b.orientation = dorientation * b.orientation

  # Compute new energy
  e_new = external_enery(b, radius)
  
  # Do Metropolies step
  mcmc = exp(-(e_new - e_current) / kT)
  rn = rand(rng, Float64)
  if rn < mcmc
    b.q_old = copy(b.q)
    b.orientation_old = deepcopy(b.orientation)   
  else
    b.q = copy(b.q_old)
    b.orientation = deepcopy(b.orientation_old) 
  end  
end


# Save last config
if (n_steps % n_save == 0) 
  write(f_config, string(1), "\n")
  println(f_config, join(b.q, " "), " ", b.orientation)
end


