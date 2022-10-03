# See PS3_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

# NOTE: If you get an error complaining about redefining a constant, you may need to kill and restart the Julia terminal
#cd("/Users/dalya/Documents/GitHub/compdjj/PS3")

using Parameters
using DelimitedFiles
using Plots

# Input - exogenous selections, also including grid parameters
@with_kw struct Input
    # Production parameters
    δ::Float64 = 0.06 # depreciation rate
    α::Float64 = 0.36 # capital share
    
    # HH life cycle parameters
    N::Int64 = 66 # number of periods agents live for
    n::Float64 = 0.011 # population growth rate
    Jr::Int64 = 46 # retirement age
    eta = open(readdlm,"ef.txt")
    
    # HH utility function parameters
    γ::Float64 = 0.42 # weight on consumption
    σ::Float64 = 2.0 # coefficient of relative risk aversion
    β::Float64 = 0.97 # discount factor
    
    # HH productivity Markov chain
    z::Array{Float64,1} = [3.0, 0.5] # high/low idiosyncratic productivity
    z_length::Int64 = length(z)
    z_dist::Array{Float64,1} = vcat(fill(3.0,2037), fill(0.5,7963))
    Π::Array{Float64,2} = [0.9261 0.07389; 0.0189 0.9811] # random productivity persistence probabilities
    Π_ergodic::Array{Float64} = [0.2037, 0.7963] # stationary distribution of the Markov chain z
    
    # Asset grid
    A::Array{Float64,1} = [0.0, 100.0] # space of asset holdings -- can't find this in the PS?
    a_length::Int64 = 200 # asset grid length, count
    a_grid::Array{Float64,1} = range(start=A[1],length=a_length,stop=A[2]) # asset grid
    
    # Prices faced by HHs
    θ::Float64 = 0.11 # proportional labor income tax
    w::Float64 = 1.05 # wage
    r::Float64 = 0.05 # interest rate
    b::Float64 = 0.2 # social security benefit
end

# Output
mutable struct Output
    valfunc::Array{Float64,3} # value function, assets x age
    polfunc::Array{Float64,3} # policy function (capital/savings)
    labfunc::Array{Float64,3} # policy function (labor choice)
    F::Array{Float64,3} # steady-state distribution of agents over age, productivity, assets
    μ_age::Array{Float64} # relative size of age cohorts
    K_SS::Float64 # steady-state capital
    L_SS::Float64 # steady-state labor
    r_SS::Float64 # steady-state interest rate
    w_SS::Float64 # steady-state wage
    b_SS::Float64 # steady-state pension benefit
end

# Struct to store transition path objects
mutable struct Transition
    T::Int64 # transition time
    
    polfunc_t::Array{Float64,4} # policy function: assets x age x state, for each time t
    labfunc_t::Array{Float64,4} # policy function: assets x age x state, for each time t
    valfunc_0::Array{Float64,3} # value function at t=0
    
    F_t::Array{Float64,4} # distribution of agents: assets x age x state, for each time t
    
    K_path::Array{Float64,1} # path of aggregate capital
    L_path::Array{Float64,1} # path of aggregate labor
    
    w_path::Array{Float64,1} # path of wages
    r_path::Array{Float64,1} # path of interest rates
    b_path::Array{Float64,1} # path of pension benefits
end

# Initialize structure to hold transition path objects
function initialize_transition(input::Input, T::Int64)
    @unpack a_length, N, Jr, z_length = input

    T = T

    polfunc_t = zeros(a_length, N, z_length, T)
    labfunc_t = zeros(a_length, N, z_length, T)
    valfunc_0 = fill(-Inf, a_length, N, z_length)

    F_t = zeros(a_length, N, z_length, T)

    K_path = zeros(T)
    L_path = zeros(T)

    w_path = zeros(T)
    r_path = zeros(T)
    b_path = zeros(T)

    return Transition(T, polfunc_t, labfunc_t, valfunc_0, F_t, K_path, L_path, w_path, r_path, b_path)
end

# Initialize structures
function Initialize(input::Input)
    @unpack a_length, N, Jr = input
    valfunc = zeros(a_length,N,2)
    polfunc = zeros(a_length,N,2)
    labfunc = zeros(a_length,Jr-1,2)
    
    # distribution over assets, age, states
    F = zeros(a_length, N, 2)

    # relative size of age cohorts
    μ_age = zeros(N)
    
    # initial guesses for steady-state capital, labor, and prices/benefits
    K_SS = 1
    L_SS = 1
    r_SS = input.r
    w_SS = input.w
    b_SS = input.b
    
    return Output(valfunc, polfunc, labfunc, F, μ_age, K_SS, L_SS, w_SS, r_SS, b_SS)
end

# Worker utility function
function uw(c::Float64, l::Float64, γ::Float64, σ::Float64)
    if c>0
        return (((c^γ)*(1-l)^(1-γ))^(1-σ))/(1-σ)
    else
        return -Inf
    end
end

# Retiree utility function
function ur(c::Float64, γ::Float64, σ::Float64)
    if c>0
        return (c^((1-σ)*γ))/(1-σ)
    else
        return -Inf
    end
end    

# Solve retiree's problem with backward induction
function retiree(input::Input, output::Output; updating::Bool=false)
    @unpack γ, σ, a_length, a_grid, N, Jr, r, b, β = input

    if updating == true
        r = output.r_SS
    end

    for j = N:-1:Jr
        for i = 1:a_length # iterate over possible asset states
            max_val = -Inf
            if j == N # last period, no real choice to make
                max_val = ur((1+r) * a_grid[i] + b, γ, σ)
            else
                for i_next = 1:a_length # iterate over assets tomorrow choices
                    v_next = output.valfunc[i_next, j+1, 1]
                    val = ur((1+r) * a_grid[i] + b - a_grid[i_next], γ, σ) + β * v_next
                    if val > max_val
                        max_val = val
                        output.polfunc[i, j, 1] = a_grid[i_next]
                        output.polfunc[i, j, 2] = a_grid[i_next]
                    end
                end
            end
            output.valfunc[i, j, 1] = max_val
            output.valfunc[i, j, 2] = max_val
        end
    end
end

# Solve worker's problem with backward induction
function worker(input::Input, output::Output; updating::Bool=false)
    @unpack γ, σ, a_length, a_grid, N, Jr, r, b, θ, β, w, eta, z, Π = input

    if updating == true
        r = output.r_SS
        w = output.w_SS
    end

    for j = (Jr-1):-1:1 # looping backwards from retirement to work
        for i_z = 1:2
            for i = 1:a_length # first loop over asset states
                max_val = -Inf # set our current max val to a low number
                for i_next = 1:a_length # loop over assets tomorrow
                    v_next = Π[i_z, 1] * output.valfunc[i_next, j+1, 1] + Π[i_z, 2] * output.valfunc[i_next, j+1, 2]
                    e = eta[j] * z[i_z]  # add age efficiency
                    l = ((γ * (1 - θ) * e * w) - ((1-γ) * ((1+r) * a_grid[i] - a_grid[i_next]))) / ((1 - θ) * w * e) # labor
                    l = min(1.0, max(0.0, l)) # make sure labor is between 0 and 1
                    val = uw(w * (1 - θ) * e * l + (1+r) * a_grid[i] - a_grid[i_next], l, γ, σ) + β * v_next
                    if val > max_val
                        max_val = val
                        output.polfunc[i, j, i_z] = a_grid[i_next]
                        output.labfunc[i, j, i_z] = l
                    end
                end
                output.valfunc[i,j,i_z] = max_val
            end
        end
    end    
end

# Compute the steady-state distribution of agents over age (j), productivity (z), asset holdings (a): F_j(z, a)
function distribution(input::Input, output::Output)
    @unpack N, Π, a_length, a_grid, n, z_length, Π_ergodic = input
    @unpack polfunc = output

    # Step 1: find relative sizes of each cohort of each age j
    μ_age = ones(N)
    for j=2:N
        μ_age[j] = μ_age[j-1] / (1+n)
    end
    μ_age = μ_age / sum(μ_age)
    output.μ_age = μ_age

    # Step 2: compute wealth distribution for each age cohort, using policy rules + distribution of prev cohort
    # age j=1: 
    @assert a_grid[1] == 0.0
    F = zeros(a_length, N, z_length)
    F[1, 1, 1] = μ_age[1] * Π_ergodic[1]
    F[1, 1, 2] = μ_age[1] * Π_ergodic[2]

    # age j>1: use decision rules + prev cohort's distribution of assets
    @time for j=2:N, i=1:a_length, s=1:z_length
        # Mathematically, this computes (for each age j, asset level (today) i, state (today) s):
        # μ_age[j] * { sum over i', s':  F[i', j-1, s'] * Indicator(polfunc[i', s', j-1] == a_grid[i]) * Π[s', s]  }
        
        temp = F[:, j - 1, :] .* (polfunc[:, j - 1, :] .== a_grid[i])
        #     # 1000x2 matrix
        #     # Each entry (i',s') is the mass of people at age (j-1) who had assets a_grid[i'] & were at state z[s'] last period,
        #     # and optimally chose a_grid[i] for the next period
        #     # i.e., temp[i', s'] = F[i', j-1, s'] * Indicator(polfunc[i', s', j-1] == a_grid[i])

        F[i, j, s] = 1/(1+n) * sum(temp' .* Π[:, s])
        #     # Multiply each term in the prev matrix by the fraction who transition from state s' (at age j-1) to state s (age j)
        #     # and then sum over s'. Finally, multiply this by the size of age cohort j
    end

    output.F = F
end

# Update current guess for steady-state aggregate labor and capital K_SS and L_SS
function K_L_update(input::Input, output::Output)
    @unpack θ, Jr, N, a_grid, a_length, z_length, eta, z = input
    #@unpack w_SS, L_SS, μ_age, F, labfunc = output

    # compute new guess for steady-state capital -- SHOULD VECTORIZE THIS
    K_next = sum(output.F .* a_grid)

    # e[z, j]: earnings in state z at age j=1, ..., (Jr - 1) (2 x 45 matrix)
    e = eta' .* z

    # compute new guess for steady-state labor
    sum1 = sum(output.labfunc[:, :, 1] .* e[1, :]' .* output.F[:, 1:(Jr-1), 1]) # dot product for state z=1
    sum2 = sum(output.labfunc[:, :, 2] .* e[2, :]' .* output.F[:, 1:(Jr-1), 2]) # dot product for state z=2
    L_next = sum1 + sum2

    return K_next, L_next
end

# Solve for a SS equilibrium:
# (1) Solve the retiree & worker problems, given current guesses of capital & labor (& thus prices and pensions).
# (2) Then compute the implied SS distribution of assets and employment states across age cohorts.
# (3) Then compute the implied aggregate stocks of capital and labor. If close to the initial guesses, stop. Otherwise,
# update guesses for K, L, prices, and pension benefits and repeat.
function solve_steady_state(input::Input, output::Output; update_factor::Float64=0.01, tol::Float64=1e-3, maxiter::Int64=10000)
    @unpack α, Jr, N, θ = input
    
    err = 1000.0

    counter = 0

    while err > tol && counter < maxiter
        println("***** Iteration $counter:")
        
        # (1) Solve retiree and worker problems with current prices
        println("Solving retiree problem")
        @elapsed retiree(input, output; updating=true)
        
        println("Solving worker problem")
        @elapsed worker(input, output; updating=true)
        
        # (2) Compute distribution of assets and productivity states across ages
        println("Solving for stationary distribution F")
        @elapsed distribution(input, output)
        
        # (3) Compute aggregate K and L implied by results from (1) and (2)
        println("Updating capital and labor guesses")
        @elapsed K_next, L_next = K_L_update(input, output)

        # Calculate difference b/w prev and current guesses for SS K and L
        err = max(abs(K_next - output.K_SS), abs(L_next - output.L_SS))
        println("Error at iteration $counter: " * string(err))
        
        # Update guesses for SS K and L
        output.K_SS = update_factor * K_next + (1 - update_factor) * output.K_SS
        output.L_SS = update_factor * L_next + (1 - update_factor) * output.L_SS

        # Update prices based on updated K and L guesses
        output.w_SS = (1-α) * (output.K_SS^α) * (output.L_SS^(-α)) # wage = MPL
        output.r_SS = α * (output.K_SS^(α - 1)) * (output.L_SS^(1 - α)) # interest = MPK
        output.b_SS = (θ * output.w_SS * output.L_SS) / sum(output.μ_age[Jr:N])

        counter += 1
    end

    if counter == maxiter
        println("Reached max number of iterations ($maxiter) for K, L convergence.")
    else
        println("Guesses for K_SS and L_SS converged in $counter iterations!")
    end
end

#---------------------------------------------------
#----Compute initial guess of the transition path, given a guess for time T
# Input: T = time at which the economy has finished transitioning to the no-Soc Sec steady state
# Output: linear paths of K and L, and the implied paths of w, r, and b.
# (Note: we assume b = the "with-Social Security" steady state value at t=1; and b=0 for all t>1.)
function K_L_path_initial(input_no_socsec::Input, output_transition::Transition)
    @unpack γ, σ, α, a_length, a_grid, N, Jr, r, b, β, z, Π, z_length, eta, θ = input_no_socsec
    @unpack T = output_transition

    # initial guesses for path of capital and labor
    output_transition.K_path = zeros(T)
    output_transition.L_path = zeros(T)
    for t=1:T
        output_transition.K_path[t] = output_w_socsec.K_SS + (t - 1) * (output_no_socsec.K_SS - output_w_socsec.K_SS) / (T - 1)
        output_transition.L_path[t] = output_w_socsec.L_SS + (t - 1) * (output_no_socsec.L_SS - output_w_socsec.L_SS) / (T - 1)
    end

    # path of prices: plug K and L paths into marginal product formulas
    output_transition.w_path = (1 - α) .* (output_transition.K_path .^ α) .* (output_transition.L_path .^ (-α)) # wage = MPL
    output_transition.r_path = α .* (output_transition.K_path .^ (α - 1)) .* (output_transition.L_path .^ (1 - α)) # interest = MPK

    # path of pension benefits: starts at initial SS, then goes to 0 in next period
    output_transition.b_path = zeros(T)
    output_transition.b_path[1] = output_w_socsec.b_SS
end

#----- Solve the HH problem in every time period, backwards from t=T
function HH_path(input_no_socsec::Input, output_no_socsec::Output, output_transition::Transition)
    @unpack γ, σ, α, a_length, a_grid, N, Jr, r, b, β, z, Π, z_length, eta, θ = input_no_socsec
    @unpack T = output_transition

    # Val func & pol funcs for t=T
    valfunc_t_next = output_no_socsec.valfunc
    output_transition.polfunc_t[:, :, :, T] = output_no_socsec.polfunc
    output_transition.labfunc_t[:, 1:(Jr-1), :, T] = output_no_socsec.labfunc

    # Solve for val func and pol funcs for t=T-1, t=T-2, ..., t=1
    @unpack r_path, w_path, b_path = output_transition

    for t=(T - 1):-1:1
        valfunc_t = fill(-Inf, a_length, N, z_length)
        # polfunc_t = zeros(a_length, N, z_length)
        # labfunc_t = zeros(a_length, Jr-1, z_length)

        # retirees
        for j=N:-1:Jr, i=1:a_length
            if j==N # oldest age cohort: don't save bc will die next period
                c = max.(0.0, (1 + r_path[t]) * a_grid[i] + b_path[t])
                valfunc_t[i, j, :] .= ur(c, γ, σ)
            else # retired cohorts who won't die next period
                # all possible continuation values, for each a'
                v_next = valfunc_t_next[:, j + 1, 1]
                
                # all possible consumption levels, for each a'
                c = max.(0.0, (1 + r_path[t]) * a_grid[i] + b_path[t] .- a_grid)
                
                # all possible lifetime values, for each a'
                vals = ur.(c, γ, σ) .+ β * v_next
                
                # Find a' that yields the max lifetime value
                max_index1 = findmax(vals)[2]
                max_index2 = findmax(vals)[2]

                # store the associated lifetime value + policy function
                # (Same value for each productivity state b/c retirees aren'
                # affected by productivity shocks z)
                valfunc_t[i, j, :] .= vals[max_index1]
                output_transition.polfunc_t[i, j, :, t] .= a_grid[max_index1]
            end
        end

        # workers
        for j=(Jr-1):-1:1, i=1:a_length, i_z=1:z_length
            # all possible expected continuation values, for each a'
            v_next = Π[i_z, 1] .* valfunc_t_next[:, j + 1, 1] + Π[i_z, 2] .* valfunc_t_next[:, j + 1, 2]

            # efficient hours worked per unit of labor
            e = eta[j] * z[i_z]  # add age efficiency

            # all possible levels of labor supply, for each a'
            l = ((γ * (1 - θ) * e * w_path[t]) .- ((1 - γ) * ((1 + r_path[t]) .* a_grid[i] .- a_grid))) ./ ((1 - θ) * w_path[t] * e) # labor
            l = min.(1.0, max.(0.0, l)) # make sure labor is between 0 and 1

            # all possible levels of consumption, for each a'
            c = max.(0.0, w_path[t] .* (1 - θ) .* e .* l .+ (1 + r_path[t]) .* a_grid[i] .- a_grid)

            # all possible lifetime values, given different choices of a' (and therefore c and l)
            vals = uw.(c, l, γ, σ) .+ β * v_next

            # find the a' which yields the max lifetime value
            # then store the associated value, savings a', and labor choice l
            max_index = findmax(vals)[2]
            valfunc_t[i, j, i_z] = vals[max_index]
            output_transition.polfunc_t[i, j, i_z, t] = a_grid[max_index]
            output_transition.labfunc_t[i, j, i_z, t] = l[max_index]
        end

        # save the results for each t; to be used on the next (t-1) iteration
        valfunc_t_next = valfunc_t
    end # end loop over t

    output_transition.valfunc_0 = valfunc_t_next
end

#-----Compute the distribution of assets across age cohorts and states, *for every t*
function distribution_path()
end

#-----Compute the path of K and L implied by the results from the above two functions
function K_L_path_update()
end

#-----Iterate the above steps until the guesses for the K,L paths converge
function K_L_path_iterate()
end