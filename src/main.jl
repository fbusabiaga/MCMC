"Main code."

include("read_input/read_input.jl")
using .read_input



print("Start\n")

# Read input
options = read_input.read_input_file("data.main")

print("options = ", options, "\n")

# Create bodies
q = zeros(Float64, 3)
print("q = ", q, "\n")

# Loop steps

# Compute energy

# Metropolies

# Save config

# Close


print("End\n")

