# See PS2_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

# NOTE: If you get an error complaining about redefining a constant, you may need to kill and restart the Julia terminal

## 0 - Libraries ##
using Parameters

## 1 - Declaring Variables ##

# Input - exogenous selections, also including grid parameters
@with_kw struct Input
    beta::Float64 = 0.9932 # discount factor
    alpha::Float64 = 1.5 # coefficient of relative risk aversion
    S::Array{Float64,1} = [1.0,0.5] # set of possible earnings, e=1,u=0.5
    P::Array{Float64,2} = [[0.97,0.5][0.03,0.5]] # transition probabilities matrix
    A::Array{Float64,1} = [-2.0,5.0] # space of asset holdings

    a_length::Int64 = 1000 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
end

# Output
mutable struct Output
    valfunc::Array{Float64,2} # value function
    polfunc::Array{Float64,2} # policy function
    q::Float64 # price of discount bonds
end

## 2 - Value Function Iteration ##

function Bellman(input::Input,output::Output)
# This function is to be called each value function iteration.
# It returns the next guess of the value function.
    @unpack valfunc = output # get last value function guess, as given
    @unpack beta, alpha, S, P, A = input # get parameters
    v_next = zeros(a_length,2) # Initializing the next guess as zeroes, 2 for number of emp. states
    # [add in for loop over emp states]
        # [add in for loop over asset grid]
            # [add in for loop over choices of next a]
    return v_next
end
