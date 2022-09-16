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
    Π::Array{Float64,2} = [0.97 0.5; 0.03 0.5] # same thing: transition probabilities matrix but to avoid issue
    A::Array{Float64,1} = [-2.0,5.0] # space of asset holdings
    price::Float64 = 0.3 # we should have a price ??
    a_length::Int64 = 1000 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
end

# Output
mutable struct Output
    valfunc::Array{Float64,3} # value function
    polfunc::Array{Float64,3} # policy function
    q::Float64 # price of discount bonds # are we choosing price or given ??
end

## 2 - Value Function Iteration ##

function Bellman(input::Input,output::Output)
# This function is to be called each value function iteration.
# It returns the next guess of the value function.
    @unpack valfunc = output # get last value function guess, as given
    @unpack beta, alpha, S, P, Π, A, a_length, a_grid = input # get parameters
    v_next = zeros(a_length,2) # Initializing the next guess as zeroes, 2 for number of emp. states
    # [add in for loop over emp states]
    for e_index = 1:2 #iterating over values of emp. states
        e = S[e_index] # which emp. state
        candidate_max = 0 # initial candidate
        # [add in for loop over asset grid]
        for a_index = 1:a_length  
            a = a_grid[a_index] 
            # [add in for loop over choices of next a]
            for ap_index = 1:a_length
                ap = a_grid[ap_index]
                c = e + a - ap*price
                if c>0
                    val = ((c^(1-alpha)-1)/(1-alpha)) + β*(Π[e_index, 1]*val_func[a_index, 1] + Π[e_index, 2]*val_func[a_index, 2])
            end
        end
    end
    return v_next
end
