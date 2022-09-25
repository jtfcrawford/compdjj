# See PS3_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

# NOTE: If you get an error complaining about redefining a constant, you may need to kill and restart the Julia terminal
#cd("/Users/dalya/Documents/GitHub/compdjj/PS3")

using Parameters
using DelimitedFiles
using Plots

# Input - exogenous selections, also including grid parameters
@with_kw struct Input
    initial_a::Float64 = 0
    δ::Float64 = 0.06 # depreciation rate
    α::Float64 = 0.36 # capital share
    N::Int64 = 66 # number of periods agents live for
    n::Float64 = 0.011 # population growth rate
    Jr::Int64 = 46 # retirement age
    θ::Float64 = 0.11 # proportional labor income tax
    γ::Float64 = 0.42 # weight on consumption
    σ::Float64 = -2.0 # coefficient of relative risk aversion
    z::Array{Float64,1} = [3, 0.5] # high/low idiosyncratic productivity
    #z_odd::Array{Float64,1} = [0.2037, 0.7963] # chances of drawing a high/low z
    β::Float64 = 0.97 # discount factor
    w::Float64 = 1.05 # wage
    r::Float64 = 0.05 # interest rate
    b::Float64 = 0.2 # social security benefit
    Π::Array{Float64,2} = [0.9261 0.07389; 0.0189 0.9811] # random productivity persistence probabilities
    #eta = int(open(readdlm,"ef.txt"))
    A::Array{Float64,1} = [-10.0, 100.0] # space of asset holdings -- can't find this in the PS?
    a_length::Int64 = 1000 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
end

# Output
mutable struct Output
    valfunc::Array{Float64,3} # value function, assets x age
    polfunc::Array{Float64,3} # policy function (capital/savings)
    labfunc::Array{Float64,3} # policy function (labor choice)
end

# Initialize structures
function Initialize(input::Input)
    @unpack a_length, N, Jr = input
    valfunc = zeros(a_length,N,2)
    polfunc = zeros(a_length,N,2)
    labfunc = zeros(a_length,Jr-1,2)
    return Output(valfunc,polfunc,labfunc)
end

# Retiree utility function
function ur(c::Float64,γ::Float64,σ::Float64)
    if c>0
        return (c^((1-σ)*γ))/(1-σ)
    else
        return -Inf
    end
end 


# Solve retiree problem (backwards)
function retiree(input::Input,output::Output)
    @unpack γ, σ, a_length, a_grid, N, Jr, r, b, β = input
    
    for j = N:-1:Jr
        for i = 1:a_length # iterate over possible asset states
            max_val = -Inf
            if j == N # last period, no real choice to make
                output.polfunc[i,j,1] = 0
                output.valfunc[i,j,1] = 0
            else
                for i_next = 1:a_length # iterate over assets tomorrow choices
                    v_next = output.valfunc[i_next,j+1,1]
                    a_next = output.polfunc[i_next,j+1,1]
                    val = ur((1+r)*a_grid[i_next]+b-a_next,γ,σ) + β*v_next
                    print(val)
                    break
                    if abs.(val - output.valfunc[i_next,j,1]) != Inf
                        if abs.(val - output.valfunc[i_next,j,1]) < 0.001
                            max_val = val
                            output.polfunc[i,j,1] = a_grid[i]
                            output.polfunc[i,j,2] = a_grid[i]
                        end
                    end
                end
            end
        end
    end    
end
#=
error = 1000.0 # Some large number so the loop doesn't immediately stop
counter = 0
while error > tol
    v_next = Bellman(input, temp_results, q)
    error = maximum(abs.(v_next .- temp_results.valfunc))
    #println(string(maximum(temp_results.valfunc)) * '\n')
    temp_results.valfunc = v_next
    counter = counter + 1
end
println("Value function converged in $counter iterations. Error: $error; Tolerance: $tol.")
# Nothing is returned -- the "result" of running this is an updated polfunc
end
=#