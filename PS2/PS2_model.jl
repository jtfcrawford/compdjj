# See PS2_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

# NOTE: If you get an error complaining about redefining a constant, you may need to kill and restart the Julia terminal

## 0 - Libraries ##
using Parameters

## 1 - Declaring Variables ##

# Input - exogenous selections, also including grid parameters
@with_kw struct Input
    β::Float64 = 0.9932 # discount factor
    α::Float64 = 1.5 # coefficient of relative risk aversion
    
    S::Array{Float64,1} = [1.0, 0.5] # set of possible earnings, e=1,u=0.5
    Π::Array{Float64,2} = [0.97 0.03; 0.5 0.5] # same thing: transition probabilities matrix but to avoid issue
    S_length::Int64 = length(S)

    A::Array{Float64,1} = [-2.0, 5.0] # space of asset holdings
    a_length::Int64 = 1000 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
    
    q_initial::Float64 = 0.993 # initial guess for price of discount bonds
    q_min::Float64 = 0 # range of possible bond prices
    q_max::Float64 = 1
end

# Output
mutable struct Output
    valfunc::Array{Float64,2} # value function
    polfunc::Array{Float64,2} # policy function
    μ::Array{Float64,2} # bond holdings / wealth distribution
    q::Float64 # equilibrium price of discount bond
end

## Initialize structures
function Initialize()
    # initialize struct with model primitives
    input = Input()

    # initialize struct with results
    valfunc = zeros(input.a_length, input.S_length)
    polfunc = zeros(input.a_length, input.S_length)

    # initial guess for stationary wealth distribution
    # try a uniform distribution (PS2 suggests something else -- maybe replace later)
    μ = fill((input.a_length * input.S_length)^(-1), input.a_length, input.S_length)

    # initial guess for market-clearing price
    q = input.q_initial

    # initialize struct with equilibrum objects to solve for
    output = Output(valfunc, polfunc, μ, q)

    return input, output
end

## 2 - Value Function Iteration ##

function Bellman(input::Input,output::Output)
# This function is to be called each value function iteration.
# It returns the next guess of the value function.
    @unpack valfunc, q = output # get last guesses for value function and bond price
    @unpack β, α, S, S_length, Π, A, a_length, a_grid = input # get parameters
    
    v_next = zeros(a_length, 2) # Initializing the next guess as zeroes, 2 for number of emp. states
    
    # for loop over emp states
    for e_index = 1:2, a_index = 1:a_length # loop over (a,s)
        
        e = S[e_index] # which emp. state
        a = a_grid[a_index]
        
        # for loop over choices of next a'
        candidate_max = -Inf # initial candidate
        for ap_index = 1:a_length
            ap = a_grid[ap_index]
            c = e + a - ap*q # consumption
            if c > 0
                val = ((c^(1 - α) - 1)/(1 - α)) + β*(Π[e_index, 1] * valfunc[ap_index, 1] + Π[e_index, 2] * valfunc[ap_index, 2])
                if val > candidate_max
                    candidate_max = val
                    output.polfunc[a_index,e_index] = ap # Updating decision rule
                end
            end
        end
        v_next[a_index,e_index] = candidate_max
    end # end loop over (a,s)
    return v_next
end

function Value_Iteration(input::Input,output::Output; tol::Float64 = 1e-5)
# This function re-runs the Bellman function until sufficient value function convergence
    error = 1000.0 # Some large number so the loop doesn't immediately stop
    counter = 0
    while error > tol
        v_next = Bellman(input,output)  
        error = maximum(abs.(v_next .- output.valfunc))
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

function T_star(input::Input,output::Output)
    @unpack S, Π, S_length, A, a_length, a_grid = input
    @unpack polfunc, μ = output

    # initial guess for next iteration of μ
    μ_next = zeros(a_length, S_length)
    
    # μ_next[i_next, j_next] = Pr(a' = a_grid[i_next], s' = S[j_next])
    # Using a law of motion, this is:
    #   whether a HH will have a'=a[i_next] tomorrow, given they today have a=a_grid[i_today] and are in state S[j_today], times
    #   the probability that a HH is today in state S[j_today] and goes to S[j_next] tomorrow, times
    #   the measure of HHs who today have a=a_grid[i_today] and are in state S[j_today]
    # ... summed over all assets a_grid[i_today] and all states S[j_today]
    for i_next=1:a_length, j_next=1:S_length
        for i_today=1:a_length, j_today=1:S_length
            μ_next[i_next,j_next] += (polfunc[i_today, j_today] == a_grid[i_next]) * Π[j_today, j_next] * μ[i_today, j_today]
        end
    end
    
    # check that μ_next approximately sums to 1 (code breaks if not)
    # sum_check = sum(μ_next)
    # @assert abs(sum(μ_next) - 1) < 0.01 "Error: sum of μ_next is $sum_check; not within 0.01 of 1."

    # update μ
    return μ_next
end

# Then we need another function that basically is the same as the value iteration function
# except we want the wealth distribution to converge
function Wealth_Dist_Iteration(input::Input, output::Output; tol::Float64=1e-5, err::Float64 = 100.0, maxiter::Int64 = 5000)
    counter = 0
    while err > tol && counter <= maxiter
        μ_next = T_star(input, output) # apply T* operator to current guess for μ*
        μ_next_sum = sum(μ_next)
        err = maximum(abs.((μ_next .- output.μ) ./ output.μ)) # update error
        println("Iteration $counter: Sum of μ_next is $μ_next_sum, err is $err")
        output.μ = μ_next
        counter += 1
    end
    if counter > maxiter
        println("Max number of iterations ($maxiter) reached for wealth distribution. Error: $err; tolerance: $tol.")
    else
        println("Wealth distribution converged in $counter iterations. Error: $err; tolerance: $tol")
    end
end

# Once we have the above set up, we can use the below price updating function and continually
# re-run value iteration and wealth distribution iteration until we find the "right" price


function adjust_price(input::Input,output::Output,tol::Float64 = 1e-5)
   excess_supply = sum(output.μ .* [input.a_grid input.a_grid])
   if abs(excess_supply) < tol

   # I think I was close on how to do this bisection but the below isn't quite right

   #elseif excess_supply > tol
   # output.q = (output.q + q_max) / 2
   #elseif excess_supply < -tol
   # output.q = (output.q + q_min) / 2
   end
end