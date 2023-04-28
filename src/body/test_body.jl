"Test the quaternion implementation."

print("Start\n")

include("body.jl")
import .body_module

# Create body
q = [0, 0]
orientation = 0
N_blobs = 4
theta = collect(range(0,2*pi,length=N_blobs+1)[1:end-1])
r_vectors = Array{Float64,2}(undef, N_blobs, 2)
r_vectors[:,1] = cos.(theta)
r_vectors[:,2] = sin.(theta)
b = body_module.body(q, orientation, r_vectors)
print("b.q           = ", b.q, "\n")
print("b.orientation = ", b.orientation, "\n")
print("b.r_vectors     = \n")
display(b.r_vectors)
print("\n")

# Update body
b.orientation = pi / 4
r_vectors = body_module.get_r_vectors(b)
print("r_vectors = \n")
display(r_vectors)
print("\n")


print("End\n")
