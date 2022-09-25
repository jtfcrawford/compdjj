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
    σ::Float64 = 2.0 # coefficient of relative risk aversion
    z::Array{Float64,1} = [3.0, 0.5] # high/low idiosyncratic productivity
    z_dist::Array{Float64,1} = vcat(fill(3.0,2037), fill(0.5,7963))
    β::Float64 = 0.97 # discount factor
    w::Float64 = 1.05 # wage
    r::Float64 = 0.05 # interest rate
    b::Float64 = 0.2 # social security benefit
    Π::Array{Float64,2} = [0.9261 0.07389; 0.0189 0.9811] # random productivity persistence probabilities
    eta = open(readdlm,"ef.txt")
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
    @unpack γ, σ, a_length, a_grid, N, Jr, r, b, β = input

    for j = N:-1:Jr
        for i = 1:a_length # iterate over possible asset states
            max_val = -Inf
            if j == N # last period, no real choice to make
                max_val = ur((1+r)*a_grid[i]+b,γ,σ)
            else
                for i_next = 1:a_length # iterate over assets tomorrow choices
                    v_next = output.valfunc[i_next,j+1,1]
                    val = ur((1+r)*a_grid[i]+b-a_grid[i_next],γ,σ) + β*v_next
                    if val > max_val
                        max_val = val
                        output.polfunc[i,j,1] = a_grid[i_next]
                        output.polfunc[i,j,2] = a_grid[i_next]
                    end
                end
            end
            output.valfunc[i,j,1] = max_val
            output.valfunc[i,j,2] = max_val
        end
    end    
end

# Solve workers problem 
function worker(input::Input,output::Output)
    @unpack γ, σ, a_length, a_grid, N, Jr, r, b, β, w, θ, eta, z_dist, Π = input
    @unpack valfunc = output
    for j = Jr:-1:1 # looping backwards from retirement to work
        for i = 1:a_length # first loop over asset states
            max_val = -Inf # set our current max val to a low number
            if j == Jr # if we are retiring, max val is the last result from retiree's problem
                max_val = output.valfunc[i,j,1]
            else
                for i_next = 1:a_length # loop over assets tomorrow
                    v_next = output.valfunc[i_next,j+1,1] # we know tomorrow's value would be this?
                    z_next = rand(z_dist,1) # update productivity
                    e = eta[j] * z_next[1]  # add age efficiency
                    l = ( (γ*(1-θ)*e*w) - ((1-γ)*((1+r)*a_grid[i]+b-a_grid[i_next])) ) / ((1-θ)*w*e) # labor
                    if z_next == 3.0 # if high productivity now, then we use Π_{HH} and Π_{LH}
                        val = ur(w*(1-θ)*e*l*(1+r)*a_grid[i]+b-a_grid[i_next],γ,σ) + β * (Π[1,1] * valfunc[i_next, j, 1] + Π[2,1] * valfunc[i_next, j, 2])
                    else # o/w we use Π_{HL} and Π_{LL}, HOWEVER!!! we're picking next period val? still have to figure out backwards here i think
                        val = ur(w*(1-θ)*e*l*(1+r)*a_grid[i]+b-a_grid[i_next],γ,σ) + β * (Π[1,2] * valfunc[i_next, j, 1] + Π[2,2] * valfunc[i_next, j, 2])
                    end
                    if val > max_val
                        max_val = val
                        output.polfunc[i,j,1] = a_grid[i_next]
                        output.polfunc[i,j,2] = a_grid[i_next]
                    end
                end
            end
            output.valfunc[i,j,1] = max_val
            output.valfunc[i,j,2] = max_val
        end
    end    
end