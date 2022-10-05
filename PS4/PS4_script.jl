# Computational Economics - Problem Set 4
# Authors: Dalya, Jackson, Jon
# 🍬👻 October 2022 👻🍬

# This is the main file - run this one.

println("Starting...")

# You may need to manually set your file path in the Julia terminal using pwd("C:\\Example\\Filepath")
# Or you can change this line of code:
#cd("C:\\Users\\jgkro\\Documents\\GitHub\\compdjj\\PS4")
#cd("/Users/dalya/Documents/GitHub/compdjj/PS4")
cd("C:\\Users\\jaxtc\\OneDrive\\Documents\\GitHub\\compdjj\\PS4")

# Bring in model and other functions
include("PS4_model.jl")

# Solve for SS with Social Security
input_w_socsec = Input(θ=0.11)
output_w_socsec = Initialize(input_w_socsec)
@time solve_steady_state(input_w_socsec, output_w_socsec; update_factor=0.5, tol=1e-3)

# Solve for SS without Social Security
input_no_socsec = Input(θ=0.0, b=0.0)
output_no_socsec = Initialize(input_no_socsec)
@time solve_steady_state(input_no_socsec, output_no_socsec; update_factor=0.5, tol=1e-3)

# Guess: economy converges to new SS after 30 periods
output_transition = initialize_transition(input_no_socsec, 30)
K_L_path_initial(input_no_socsec, output_transition)

# Solve HH problem for each period t=T, T-1, ..., 2, 1
# Get policy functions at every t, and value function at t=1
HH_path(input_no_socsec, output_no_socsec, output_transition)

distribution_path(input_no_socsec, output_transition)

# Exercise 1-2 Plots
using Plots
plot(Array(range(1,length(output_transition.K_path))), output_transition.K_path, 
    title="K path", 
    xlabel="t", ylabel="K", legend=:topleft, linewidth=2)
savefig("PS04_Kpath.png")
plot(Array(range(1,length(output_transition.K_path))), output_transition.L_path, 
    title="L path", 
    xlabel="t", ylabel="L", legend=:topleft, linewidth=2)
savefig("PS04_Lpath.png")
plot(Array(range(1,length(output_transition.K_path))), output_transition.w_path, 
    title="w path", 
    xlabel="t", ylabel="w", legend=:topleft, linewidth=2)
savefig("PS04_wpath.png")
plot(Array(range(1,length(output_transition.K_path))), output_transition.r_path, 
    title="r path", 
    xlabel="t", ylabel="r", legend=:topleft, linewidth=2)
savefig("PS04_rpath.png")