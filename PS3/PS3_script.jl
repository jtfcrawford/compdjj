# Computational Economics - Problem Set 3
# Authors: Dalya, Jackson, Jon
# September 2022

# This is the main file - run this one.

println("Starting...")

# You may need to manually set your file path in the Julia terminal using pwd("C:\\Example\\Filepath")
# Or you can change this line of code:
#cd("C:\\Users\\jgkro\\Documents\\GitHub\\compdjj\\PS2")

# Bring in model and other functions
include("PS3_model.jl")

# Initialize input (primitives) and output (solutions) structures
input, output = Initialize()

println("Done!")