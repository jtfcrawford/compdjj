# Computational Economics - Problem Set 2
# Authors: Dalya Elmalt, Jackson Crawford, Jon Kroah
# Submitted 23 September 2022

# This is the main file - run this one.
cd("C:\\Users\\jgkro\\Documents\\GitHub\\compdjj\\PS2")

println("Starting...")

#---------------- Bring in model and other functions
include("PS2_model.jl")

#---------------- Solve the model
input, eqm = Initialize() # initialize primitives
Solve_Model(input, eqm)

#---------------- Unpack solutions for plotting
using Plots
@unpack valfunc, polfunc, μ, q = eqm
@unpack a_grid, a_length, S_length, S = input

#---------------- Report equilibrium objects
# equilibrium price / interest rate
println("Equilibrium price is $q. Implied real interest rate is $(1/q - 1).")

# conditional wealth distributions
μ_conditional = copy(μ)
μ_conditional[: , 1] = μ_conditional[:, 1] / sum(μ_conditional[:, 1])
μ_conditional[: , 2] = μ_conditional[:, 2] / sum(μ_conditional[:, 2])
plot(input.a_grid, μ_conditional, labels=["Employed" "Unemployed"], title="Distribution of Asset Holdings in Each State: μ*(a|s)", xlabel="a", ylabel="Pr(a|s)", legend=:topleft)
savefig("PS02_distributions.png")

# value functions
plot(a_grid, valfunc, labels=["Employed" "Unemployed"], title="Value Functions", xlabel="a", ylabel="V(a,s)", legend=:topleft, linewidth=2)
savefig("PS02_valfunc.png")

# policy functions
plot(a_grid, polfunc, labels=["Employed" "Unemployed"], title="Policy Functions", xlabel="a", ylabel="g(a,s)", legend=:topleft, linewidth=2)
plot!(a_grid, a_grid, color=:black, linestyle=:dash, label="g(a,s) = a")
savefig("PS02_polfunc.png")

#---------------- Lorenz Curve
cumul_pct_pop, cumul_pct_wealth, gini = Lorenz_Curve(input, eqm)

plot(cumul_pct_pop, cumul_pct_wealth, legend=:topleft, label="Lorenz Curve", title="Lorenz Curve", xlabel="Cumulative % of Population", ylabel="Cumulative % of Wealth", linewidth=2)
plot!(cumul_pct_pop, cumul_pct_pop, color=:black, linestyle=:dash, label="Perfect Equality")
plot!(cumul_pct_pop, zeros(length(cumul_pct_pop)), color=:black, label="")
savefig("PS02_lorenz.png")

println("Gini index is $gini")

#---------------- Welfare Comparisons
W_FB = Compute_Welfare_Planner(input)

W_INC = Compute_Welfare_Incomplete(eqm)

λ = Compute_Consumption_Equiv(input, eqm, W_FB)

WG = Compute_Welfare_Gains(eqm, λ)

FB_voteshare = Compute_FB_Votes(eqm, λ)

plot(a_grid, λ, labels=["Employed" "Unemployed"], title="Consumption Equivalents by Employment State", linewidth=2, xlabel="a", ylabel="Consumption Equiv. λ(a,s)")
savefig("PS02_consumption_equivs.png")

println("W_FB (first-best welfare): $W_FB")
println("W_INC (incomplete markets welfare): $W_INC")
println("WG (economywide welfare gains): $WG")
println("Fraction who would vote for complete markets: $FB_voteshare")

println("*** Done! ***")