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
    A::Float64 = 0.005
    cf::Float64 = 10.0
    cf::Float64 = 5.0
end

mutable struct Output
    p::Float64
    val_func::Array{Float64,1}
    exit_func::Array{Float64,1}
end

function initialize_output(input::Input)
    

# psθn^(θ-1)=1
# n^(θ-1)=1/(psθ)
# n=(1/psθ)^(1/(θ-1))
function labor_choice(input::Input,p::Float64,s::Float64)
    @unpack θ = input
    return max(0,(1/p*s*θ)^(1/(θ-1)))
end

function profit(input::Input,p::Float64,s::Float64,n::Float64)
    @unpack θ, cf = input
    return p*s*(n^θ)-n-p*cf
end

function exit_choice(input::Input,output::Output)
    @unpack s_array, F, β = input
    v_next = zeros(length(s_array))
    ex_next = zeros(length(s_array))
    for s_i in range(1,length(s_array))
        s = s_array[s_i]
        n = labor_choice(input,output.p,s)
        π = profit(input,output.p,s,n)
        exit_val = π
        stay_val = π + β * output.val_func * F[s_i,:]
        if exit_val > stay_val
            v_next[s_i] = exit_val
            ex_next[s_i] = 1
        else
            v_next[s_i] = stay_val
            ex_next[s_i] = 0
        end
    end
    return v_next, ex_next
end

function value_iteration(input::Input,output::Output,tol=1e-4)
    count = 0
    while max_err > tol && count < 1000
        v_next, ex_next = exit_choice(input,output)
        max_err = max(sum(abs(v_next-output.val_func)),sum(abs(ex_next-output.exit_func)))
        output.val_func = v_next
        output.exit_func = ex_next
        count = count + 1
    end
    println("Converged in $count iterations.")
end