"Test the quaternion implementation."

print("Start\n")

import Random as rnd

include("quaternion.jl")

# Create identity quaternion
print("Unit quaternion\n")
one = [1, 0, 0, 0]
q = quaternion(one)
R = rotation_matrix(q)
print("q = ", q, "\n")
print("R = \n")
display(R)
print("\n")
print("\n\n\n")


# Create quaternion for 90 degrees rotation
print("Quaternion rotated 90 degrees\n")
q = quaternion(1, [1,0,0])
R = rotation_matrix(q)
print("q = ", q, "\n")
print("R = \n")
display(R)
print("\n")
print("\n\n\n")


# Create random quaternion
print("Random quaternion\n")
rng = rnd.MersenneTwister(1234)
theta = randn(Float64, 4)
print("theta = ", theta, "\n")
q = quaternion(theta)
print("q = ", q, "\n")
print("\n\n\n")






print("End\n")
