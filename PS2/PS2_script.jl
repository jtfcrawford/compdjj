# Computational Economics - Problem Set 2
# Authors: Dalya, Jackson, Jon
# September 2022

# This is the main file - run this one.

println("Starting...")

# You may need to manually set your file path in the Julia terminal using pwd("C:\\Example\\Filepath")

# Bring in model and other functions
include("PS2_model.jl")

# Initialize
input = Input()
@unpack a_length, a_grid = Input()
valfunc = zeros(a_length,2)
polfunc = zeros(a_length,2)
mu = zeros(a_length,2)
q = 0.5
q_min = 0
q_max = 1
output = Output(valfunc,polfunc,mu)

# Solve
Value_Iteration(input,output)

# Plotting
using Plots
plot(a_grid, output.valfunc,labels=["Employed" "Unemployed"])
savefig("PS2\\valfunc.png")
plot(a_grid, output.polfunc,labels=["Employed" "Unemployed"])
savefig("PS2\\polfunc.png")

println("Done!")