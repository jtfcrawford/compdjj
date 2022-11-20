# Computational Economics - Problem Set 6

# This is the main file - run this one.

using LinearAlgebra

println("Starting...")

# You may need to manually set your file path in the Julia terminal using pwd("C:\\Example\\Filepath")
# Or you can change this line of code:
#cd("C:\\Users\\jgkro\\Documents\\GitHub\\compdjj\\PS4")
#cd("/Users/dalya/Documents/GitHub/compdjj/PS4")
cd("C:\\Users\\jaxtc\\OneDrive\\Documents\\GitHub\\compdjj\\PS7")

# Bring in model and other functions
include("PS7_functions.jl")

params = Params()

true_data = simulate_true(params,1)

#et = e_draw(1.0,10,params)

true_moments = moments_array(true_data,true)

ρ_grid = range(start=0.35,length=10,stop=0.65)

σ_grid = range(start=0.8,length=10,stop=1.2)

for ρ in ρ_grid
    for σ in σ_grid
        objective(true_moments,ρ,σ,Matrix{Float64}(I,2,2),1)
    end
end