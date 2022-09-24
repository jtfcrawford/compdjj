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
input = Input()
output = Initialize(input)

retiree(input,output)

# value functions
@unpack a_grid = input
plot(a_grid, [output.valfunc[:,50,1],output.valfunc[:,60,1]], 
    labels=["N=50" "N=60"], title="Value Functions", 
    xlabel="a", ylabel="V(a,s)", legend=:topleft, linewidth=2)
savefig("PS03_valfunc_5060.png")

# policy functions
plot(a_grid, [output.polfunc[:,50,1],output.polfunc[:,60,1]], 
    labels=["N=50" "N=60"], 
    title="Policy Functions", xlabel="a", ylabel="g(a,s)", 
    legend=:topleft, linewidth=2)
savefig("PS03_polfunc_5060.png")

println("Done!")