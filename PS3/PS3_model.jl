# See PS3_script.jl which uses this file
# This file includes the economic model itself as well as other functions.

# NOTE: If you get an error complaining about redefining a constant, you may need to kill and restart the Julia terminal
#cd("/Users/dalya/Documents/GitHub/compdjj/PS3")

using Parameters
using DelimitedFiles

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
    #zH::Float64 = 3.0 # high idiosyncratic productivity
    #zL::Float64 = 0.5 # low idiosyncratic productivity
    z::Array{Float64,1} = [3, 0.5]
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

#=
function Bellman(input::Input,output::Output)
    @unpack valfunc
    @unpack β, z, γ,σ, eta, Π, A, a_length, a_grid = input 
    v_next = zeros(a_length, 1) 
    for z_index = 1:2, a_index = 1:a_length # loop over zH and zL, and over assets
        z = zH 
        z_curr = z[z_index] # still have to figure out how people receive their initial shock
        a = a_grid[a_index]
        candidate_max = -Inf 
        for ap_index = 1:a_length
            ap = a_grid[ap_index]
            l = (γ*w*(1+θ)*z_curr*eta[j] - (1-γ)*((1+r)a)) / (w*(1-θ)*z_curr*eta[j])
            val = uw(w*(1-θ)*z_curr*eta[j]*l*a_grid[i],l,γ,σ) + β*(Π[z_index, 1] * valfunc[ap_index, 1] + Π[z_index, 2] * valfunc[ap_index, 2])
            if val > candidate_max
                candidate_max = val
                output.polfunc[a_index,e_index] = ap # Updating decision rule
            end
        end
        v_next[a_index] = candidate_max
    end 
    return v_next
end

function Value_Iteration(input::Input,output::Output; tol::Float64 = 1e-5)
    # This function re-runs the Bellman until sufficient convergence
        error = 1000.0 
        counter = 0
        while error > tol
            v_next = Bellman(input,output)  
            error = maximum(abs.(v_next .- output.valfunc))
            output.valfunc = v_next
            counter = counter + 1
        end
        print("Iterations: " * string(counter) * '\n')
        # Nothing is returned -- the "result" of running this is an updated polfunc
    end
=#
    

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

#=
# Solve worker problem (backwards, after retiree problem)
function worker(input::Input,output::Output)
    @unpack γ, σ, eta, a_length, a_grid = Input()

    for j = Jr:-1:1
        for i = 1:a_length # iterate over possible asset states
            if j == Jr # retirement time, switches to solve the retiree problem?
                retiree()
            else
                error = 1000.0 
                counter = 0
                while error > tol
                    v_next = Bellman(input,output)  
                    error = maximum(abs.(v_next .- output.valfunc))
                    output.valfunc = v_next
                    counter = counter + 1
                end
                print("Iterations: " * string(counter) * '\n')
                # Nothing is returned -- the "result" of running this is an updated polfunc
            end
        end 
    end

end
=#