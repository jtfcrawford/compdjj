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

# Solve the model
input, eqm = Initialize() # initialize primitives
Solve_Model(input, eqm)


using Plots
plot(1:length(q_guesses), q_guesses)

######Plotting the wealth distribution (TEST)
using Plots
@unpack μ = output.μ
μ_conditional = copy(μ)
μ_conditional[: , 1] = μ_conditional[:, 1] / sum(μ_conditional[:, 1])
μ_conditional[: , 2] = μ_conditional[:, 2] / sum(μ_conditional[:, 2])
plot(input.a_grid, μ_conditional, labels=["Employed" "Unemployed"], title="Wealth Distributions for Each State: μ(a|s)")

# Plotting
using Plots
plot(input.a_grid, output.valfunc,labels=["Employed" "Unemployed"])
# savefig("PS2\\valfunc.png")
savefig("valfunc.png")
plot(input.a_grid, output.polfunc,labels=["Employed" "Unemployed"])
# savefig("PS2\\polfunc.png")
savefig("polfunc.png")

println("Done!")