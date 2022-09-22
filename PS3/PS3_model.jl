# See PS3_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

# NOTE: If you get an error complaining about redefining a constant, you may need to kill and restart the Julia terminal

using Parameters

# Input - exogenous selections, also including grid parameters
@with_kw struct Input
    δ::Float64 = 0.06 # depreciation rate
    α::Float64 = 0.36 # capital share
    N::Int64 = 66 # number of periods agents live for
    n::Float64 = 0.011 # population growth rate
    Jr::Int64 = 46 # retirement age
    θ::Float64 = 0.11 # proportional labor income tax
    γ::Float64 = 0.42 # weight on consumption
    σ::Float64 = 2.0 # coefficient of relative risk aversion
    zH::Float64 = 3.0 # high idiosyncratic productivity
    zL::Float64 = 0.5 # low idiosyncratic productivity
    β::Float64 = 0.97 # discount factor

    A::Array{Float64,1} = [-10.0, 100.0] # space of asset holdings -- can't find this in the PS?
    a_length::Int64 = 1000 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
end

# Output
mutable struct Output
    valfunc::Array{Float64,2} # value function, assets x age
    polfunc::Array{Float64,2} # policy function (capital/savings)
    labfunc::Array{Float64,2} # policy function (labor choice)
end

# Initialize structures
function Initialize()
    # initialize struct with model primitives
    input = Input()

    return input
end

# Worker utility function
function uw(c::Float64,l::Float64,γ::Float64,σ::Float64)
    if c>0
        return (((c^γ)*(1-l)^(1-γ))^(1-σ))/(1-σ)
    else
        return -Inf
    end
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
    @unpack γ, σ, a_length, a_grid = Input()

    for j = N:-1:Jr
        for i = 1:a_length # iterate over possible asset states
            if j == N # last period, no real choice to make
                max_val = ur((1+r)*a_grid[i]+b,γ,σ)
            else
                for i_next = 1:a_length # iterate over assets tomorrow choices
                    v_next = output.valfunc[i_next,j+1]
                end
            end
            output.valfunc[i,j] = max_val
        end
    end
    
end

# Solve worker problem (backwards, after retiree problem)
function worker(input::Input,output::Output)
end