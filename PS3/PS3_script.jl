# Computational Economics - Problem Set 3
# Authors: Dalya, Jackson, Jon
# September 2022

# This is the main file - run this one.

println("Starting...")

# You may need to manually set your file path in the Julia terminal using pwd("C:\\Example\\Filepath")
# Or you can change this line of code:
cd("C:\\Users\\jgkro\\Documents\\GitHub\\compdjj\\PS3")
#cd("/Users/dalya/Documents/GitHub/compdjj/PS3")
# cd("C:\\Users\\jaxtc\\OneDrive\\Documents\\GitHub\\compdjj\\PS3")

# Bring in model and other functions
include("PS3_model.jl")
#include("test.jl")

# Initialize input (primitives) and output (solutions) structures
input = Input()
output = Initialize(input)

retiree(input,output)
worker(input,output)
println("Solved retiree and worker problems.")

# compute stationary distribution F_j(a,z)
distribution(input,output)

# iterate to get SS capital, labor, prices, pension benefit
K_L_iterate(input, output)

# value functions
using Plots
@unpack a_grid = input
plot(a_grid, [output.valfunc[:,20,1],output.valfunc[:,20,2],
    output.valfunc[:,50,1]], 
    labels=["N=20 (high z)" "N=20 (low z)" "N=50"], title="Value Functions", 
    xlabel="a", ylabel="V(a,s)", legend=:topleft, linewidth=2)
savefig("PS03_valfunc_2050.png")

# policy functions
plot(a_grid, [output.polfunc[:,20,1],output.polfunc[:,20,2],
    output.polfunc[:,50,1]], 
    labels=["N=20 (high z)" "N=20 (low z)" "N=50"], 
    title="Savings Functions", xlabel="a", ylabel="g(a,s)", 
    legend=:topleft, linewidth=2)
savefig("PS03_polfunc_2050.png")

println("Done!")