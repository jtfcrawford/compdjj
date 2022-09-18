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
    Π::Array{Float64,2} = [0.97 0.5; 0.03 0.5] # same thing: transition probabilities matrix but to avoid issue
    A::Array{Float64,1} = [-2.0,5.0] # space of asset holdings
    a_length::Int64 = 1000 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
    q::Float64 = 0.5 # price of discount bonds
    q_min::Float64 = 0
    q_max::Float64 = 1
end

# Output
mutable struct Output
    valfunc::Array{Float64,2} # value function
    polfunc::Array{Float64,2} # policy function
    mu::Array{Float64,2} # bond holdings / wealth distribution
end

## 2 - Value Function Iteration ##

function Bellman(input::Input,output::Output)
# This function is to be called each value function iteration.
# It returns the next guess of the value function.
    @unpack valfunc = output # get last value function guess, as given
    @unpack beta, alpha, S, Π, A, a_length, a_grid = input # get parameters
    v_next = zeros(a_length,2) # Initializing the next guess as zeroes, 2 for number of emp. states
    # for loop over emp states
    for e_index = 1:2 #iterating over values of emp. states
        e = S[e_index] # which emp. state
        candidate_max = -100000 # initial candidate
        # for loop over asset grid
        for a_index = 1:a_length  
            a = a_grid[a_index] 
            # for loop over choices of next a
            for ap_index = 1:a_length
                ap = a_grid[ap_index]
                c = e + a - ap*q # consumption
                if c>0 
                    val = ((c^(1-alpha)-1)/(1-alpha)) + beta*(Π[e_index, 1]*valfunc[a_index, 1] + Π[e_index, 2]*valfunc[a_index, 2])
                    if val > candidate_max
                        candidate_max = val
                        output.polfunc[a_index,e_index] = ap # Updating decision rule
                    end
                end
            end
            v_next[a_index,e_index] = candidate_max
        end
    end
    return v_next
end

function Value_Iteration(input::Input,output::Output,tol::Float64 = 1e-5)
# This function re-runs the Bellman function until sufficient value function convergence
    error = 1000.0 # Some large number so the loop doesn't immediately stop
    counter = 0
    while error > tol
        v_next = Bellman(input,output)  
        error = abs(maximum(v_next-output.valfunc)) / abs(maximum(v_next))
        #print(string(maximum(output.valfunc)) * '\n')
        output.valfunc = v_next
        counter = counter + 1
    end
    print("Error: " * string(error) * '\n') # print statements for current troubleshooting
    print("Iterations: " * string(counter) * '\n')
    # Nothing is returned -- the "result" of running this is an updated polfunc
end

## 3 - Bond Holding / Wealth Distribution ##

# See "Notes for PS2" on Dean's website, we want a new transition matrix that answers:
# "What is the probability that an employed person (in state e) with assets worth one
# year's income (in state a=1) is in the third wealth quintile next year?" for example.

function new_Π(input::Input,output::Output)

end

# Then we need another function that basically is the same as the value iteration function
# except we want the wealth distribution to converge

function dist_iteration(input::Input,output::Output)

end

# Once we have the above set up, we can use the below price updating function and continually
# re-run value iteration and wealth distribution iteration until we find the "right" price

function adjust_price(input::Input,output::Output,tol::Float64 = 1e-5)
   excess_supply = sum(output.mu .* [a_grid a_grid])
   if abs(excess_supply) < tol

   # I think I was close on how to do this bisection but the below isn't quite right

   #elseif excess_supply > tol
   # output.q = (output.q + q_max) / 2
   #elseif excess_supply < -tol
   # output.q = (output.q + q_min) / 2
   end
end
