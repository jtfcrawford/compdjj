# See PS2_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

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
    
    μ_initial::Array{Float64,2} = fill((a_length * S_length)^(-1), a_length, S_length) # initial guess for wealth dist -- uniform over (a,s)

    q_initial::Float64 = β # initial guess for price of discount bonds
    q_min_initial::Float64 = 0 # range of possible bond prices
    q_max_initial::Float64 = 1
end

# Struct to hold intermediate outputs: these reset every time we update q
mutable struct Temp_Results
    valfunc::Array{Float64,2} # value function
    polfunc::Array{Float64,2} # policy function
    μ::Array{Float64,2} # bond holdings / wealth distribution
end

# Struct to hold final equilibrium objects: HH val fns + policy fns; wealth dist; and market-clearing price
mutable struct Equilibrium
    valfunc::Array{Float64,2}
    polfunc::Array{Float64,2}
    μ::Array{Float64,2}
    q::Float64
    excess_demand::Float64
end

## Initialize structures
function Initialize()
    # struct with primitives/starting values
    input = Input()

    # struct with equilibrium objects (start with zeros/initial guesses)
    valfunc = zeros(input.a_length, input.S_length)
    polfunc = zeros(input.a_length, input.S_length)
    μ = input.μ_initial
    q = input.q_initial
    excess_demand = 1000.0
    eqm = Equilibrium(valfunc, polfunc, μ, q, excess_demand)

    return input, eqm
end

function Reset_Temp_Results(input::Input)
    # initialize struct with results
    valfunc = zeros(input.a_length, input.S_length)
    polfunc = zeros(input.a_length, input.S_length)

    # initial guess for stationary wealth distribution
    # try a uniform distribution (PS2 suggests something else -- maybe replace later)
    μ = input.μ_initial

    # initialize struct with equilibrum objects to solve for
    temp_results = Temp_Results(valfunc, polfunc, μ)

    return temp_results
end

## 2 - Value Function Iteration ##

# This function is to be called each value function iteration.
# In: primitives (Input), current val func / policy func / wealth dist (Temp_Results), current price (q)
# Out: the next guess of the value function.
function Bellman(input::Input, temp_results::Temp_Results, q::Float64)
    @unpack valfunc = temp_results # get last guesses for value function and bond price
    @unpack β, α, S, S_length, Π, A, a_length, a_grid = input # get parameters
    
    v_next = zeros(a_length, 2) # Initializing the next guess as zeroes, 2 for number of emp. states
    
    # for loop over emp states
    for e_index = 1:2
        ap_lower = 1 # reset starting point for a' search
        e = S[e_index] # which emp. state
        for a_index = 1:a_length # loop over (a,s)     
            a = a_grid[a_index]
            
            # for loop over choices of next a'
            candidate_max = -Inf # initial candidate
            for ap_index = ap_lower:a_length
                ap = a_grid[ap_index]
                c = e + a - ap*q # consumption
                if c > 0
                    val = ((c^(1 - α) - 1)/(1 - α)) + β*(Π[e_index, 1] * valfunc[ap_index, 1] + Π[e_index, 2] * valfunc[ap_index, 2])
                    if val > candidate_max
                        candidate_max = val
                        temp_results.polfunc[a_index, e_index] = ap # Updating decision rule
                        ap_lower = ap_index # Exploiting monotonicity of policy function
                    end
                end
            end
            v_next[a_index,e_index] = candidate_max
        end # end loop over a
    end # end loop over s
    return v_next
end

function Value_Iteration(input::Input, temp_results::Temp_Results, q::Float64; tol::Float64 = 1e-4)
# This function re-runs the Bellman function until sufficient value function convergence
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
## 3 - Bond Holding / Wealth Distribution ##

#= This function updates the guess for the stationary wealth distribution μ* (i.e., applies the T* operator).
Procedure: for all possible asset holdings tomrorow a[i_next],
and all possible emp states tomrorow S[j_next],
we want to caluclate μ_next[i_next, j_next] := Pr(a' = a_grid[i_next], s' = S[j_next]).
Using a law of motion, this is:
  whether a HH will have a'=a[i_next] tomorrow, given they today have a=a_grid[i_today] and are in state S[j_today], times
  the probability that a HH is today in state S[j_today] and goes to S[j_next] tomorrow, times
  the measure of HHs who today have a=a_grid[i_today] and are in state S[j_today]
... summed over all assets a_grid[i_today] and all states S[j_today] =#
function T_star(input::Input, temp_results::Temp_Results)
    @unpack S, Π, S_length, A, a_length, a_grid = input
    @unpack polfunc, μ = temp_results

    # fill in next guess for μ*
    μ_next = zeros(a_length, S_length)
    for i_next=1:a_length, j_next=1:S_length
        for i_today=1:a_length, j_today=1:S_length
            μ_next[i_next, j_next] += (polfunc[i_today, j_today] == a_grid[i_next]) * Π[j_today, j_next] * μ[i_today, j_today]
        end
    end

    return μ_next
end

# Then we need another function that basically is the same as the value iteration function
# except we want the wealth distribution to converge
function Wealth_Dist_Iteration(input::Input, temp_results::Temp_Results; tol::Float64=1e-5, err::Float64 = 100.0, maxiter::Int64 = 5000)
    counter = 0
    while err > tol && counter <= maxiter
        μ_next = T_star(input, temp_results) # apply T* operator to current guess for μ*
        μ_next_sum = sum(μ_next)
        err = maximum(abs.(μ_next - temp_results.μ)) # update error (sup norm)
        # if mod(counter, 100) == 0
        #     println("Wealth dist iteration $counter: Sum of μ_next is $μ_next_sum, err is $err")
        # end
        temp_results.μ = μ_next
        counter += 1
    end
    if counter > maxiter
        println("Max number of iterations ($maxiter) reached for wealth distribution. Error: $err; tolerance: $tol.")
    else
        μ_sum = sum(temp_results.μ)
        println("Wealth distribution converged in $counter iterations. Sum: $μ_sum; Error: $err; Tolerance: $tol")
    end
end

######################
# Continually re-run value iteration and wealth distribution iteration until we find the "right" price
function Solve_Model(input::Input, eqm::Equilibrium; mc_maxiter::Int64=100, mc_tol::Float64=1e-4)
    @unpack q_initial, q_min_initial, q_max_initial = input

    q = q_initial # initial guess for price
    q_min = q_min_initial
    q_max = q_max_initial
    market_cleared = 0 # set MC flag to 0
    mc_counter = 0 # count number of price guesses we'll go through

    while market_cleared == 0 && mc_counter < mc_maxiter
        println("**************************************")
        println("Guess #$mc_counter for price: q = $q")
        
        #reset value & policy functions and wealth distribution
        temp_results = Reset_Temp_Results(input)

        #given q, solve for value & policy functions
        Value_Iteration(input, temp_results, q)

        #given q and VFI results, solve for stationary wealth distribution
        Wealth_Dist_Iteration(input, temp_results)
      
        #given the above, compute excess demand and adjust price up/down.
        #excess demand for asset = sum of savings decisions weighted by μ.
        excess_demand = sum(temp_results.μ .* temp_results.polfunc)

        # store current guess for equilibrium
        eqm.valfunc = temp_results.valfunc
        eqm.polfunc = temp_results.polfunc
        eqm.μ = temp_results.μ
        eqm.q = q
        eqm.excess_demand = excess_demand

        # check if market clears; if not, raise or lower price.
        # scale the price adjustment by excess demand (smaller excess demand = closer to eqm, so adjust price less)
        if abs(excess_demand) < mc_tol
            println("Market cleared!")
            market_cleared = 1
        elseif excess_demand >= mc_tol
            q_min = q # update lower bound of q guesses to be the old q guess
            q_new = q + min(1, abs(excess_demand)) * (q_max - q) / 2 # update guess for q
            q = q_new
            mc_counter += 1
            println("Excess demand is positive ($excess_demand) and outside tolerance ($mc_tol). Raising price to $q_new")
        elseif excess_demand <= -mc_tol
            q_max = q # update upper bound of q guesses to be the old q guess
            q_new = q - min(1, abs(excess_demand)) * (q - q_min) / 2 # update guess for q
            q = q_new
            mc_counter += 1
            println("Excess demand is negative ($excess_demand) and outside tolerance ($mc_tol). Lowering price to $q_new")
        end
    end

    if mc_counter < mc_maxiter
        println("Found equilibrium price in $mc_counter guesses: q = $(eqm.q)")
        println("Excess demand = $(eqm.excess_demand); tolerance = $mc_tol.")
    else
        println("Reached max number of guesses for price ($mc_maxiter). Last guess was $(eqm.q).")
        println("Excess demand = $(eqm.excess_demand); tolerance = $mc_tol.")
    end
end
