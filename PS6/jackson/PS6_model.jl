# Functions for PS6_script.jl

@with_kw struct Input
    s_array::Array{Float64,1} = [3.98e-4 3.58 6.82 12.18 18.79]
    F::Array{Float64,2} = [0.6598 0.2600 0.0416 0.0331 0.0055;
                           0.1997 0.7201 0.0420 0.0326 0.0056;
                           0.2000 0.2000 0.5555 0.0344 0.0101;
                           0.2000 0.2000 0.2502 0.3397 0.0101;
                           0.2000 0.2000 0.2500 0.3400 0.0100]
    θ::Float64 = 0.64
    β::Float64 = 0.80
end

mutable struct Output
    p::Float64
end

# psθn^(θ-1)=1
# n^(θ-1)=1/(psθ)
# n=(1/psθ)^(1/(θ-1))
function labor_choice(input::Input,p::Float64,s::Float64)
    @unpack θ = input
    return max(0,(1/p*s*θ)^(1/(θ-1)))
end

function exit_choice(input::Input,output::Output)
    @unpack s_array, F = input
    for s_i in range(1,length(s_array))
        s = s_array[s_i]
        n = labor_choice(input,output.p,s)
    end
end