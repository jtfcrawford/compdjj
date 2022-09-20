# Computational Economics - Problem Set 2
# Authors: Dalya, Jackson, Jon
# September 2022

# This is the main file - run this one.

println("Starting...")

# You may need to manually set your file path in the Julia terminal using pwd("C:\\Example\\Filepath")
# Or you can change this line of code:
cd("C:\\Users\\jgkro\\Documents\\GitHub\\compdjj\\PS2")

# Bring in model and other functions
include("PS2_model.jl")

# Initialize input (primitives) and output (solutions) structures
input, output = Initialize()

# Find value function & policy function
Value_Iteration(input,output)

# Find stationary wealth distribution
Wealth_Dist_Iteration(input, output)

# Plotting
using Plots
plot(input.a_grid, output.valfunc,labels=["Employed" "Unemployed"])
# savefig("PS2\\valfunc.png")
savefig("valfunc.png")
plot(input.a_grid, output.polfunc,labels=["Employed" "Unemployed"])
# savefig("PS2\\polfunc.png")
savefig("polfunc.png")

println("Done!")