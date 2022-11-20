# Functions for PS7_script.jl
using Parameters
using Random
using Distributions
using DSP

@with_kw struct Params

    ρ_0::Float64 = 0.5
    σ_0::Float64 = 1.0
    μ_ϵ::Float64 = 0.0
    x_0::Float64 = 0.0
    T::Int64 = 200

end

function simulate_true(params::Params,H::Int64)
    @unpack T, x_0, ρ_0, μ_ϵ, σ_0 = params

    data = zeros(H,T)

    for i in 1:H
        for j in 1:T
            if j == 1
                data[i,j] = 0  
            else
                ϵ = rand(Normal(μ_ϵ,σ_0^2))
                data[i,j] = ρ_0*data[i,j-1] + ϵ
            end
        end
    end
    return data
end

function simulate_model(params::Params,H::Int64,σ::Float64,ρ::Float64)
    @unpack T = params

    data = zeros(H,T)

    for i in 1:H
        for j in 1:T
            if j == 1
                data[i,j] = 0  
            else
                ϵ = rand(Normal(0,σ^2))
                data[i,j] = ρ*data[i,j-1] + ϵ
            end
        end
    end
    return data
end

"""
function e_draw(var::Float64,H::Int64,params::Params)
    @unpack T = params
    et = zeros(H,T)
    for i in 1:H
        for j in 1:T
            et[i,j] = rand(Normal(0,var))
        end
    end
    return et
end
"""

function moments_array(data::Array{Float64,2},no_ac::Bool)
    H = size(data)[1]
    m_array = zeros(H,3)
    for i in 1:H
        m_array[i,1] = mean(data[i,:])
        mmean = mean(data[i,:])
        m_array[i,2] = var(data[i,:])
        autocor = zeros(size(data)[2]-1)
        for j in 2:size(data)[2]
            autocor[j-1] = (data[i,j] - mmean) * (data[i,j-1] - mmean)
        end
        m_array[i,3] = mean(autocor)
    end
    if no_ac == true
        return m_array[:,1:2]
    else
        return m_array
    end
end

function objective(true_moments::Array{Float64,2},ρ::Float64,σ::Float64,W::Array{Float64,2},H::Int64)
    model = simulate_model(params,H,σ,ρ) 
    model_moments = moments_array(model,true)

    return (true_moments - model_moments) * W * transpose(true_moments - model_moments)
end