using DIVAnd_HFRadar
using DIVAnd
using Test
using LinearAlgebra

const m = 10
const n = 11
sv = statevector((trues((m-1,n)),trues((m,n-1))))

function fun(x,δx)
    Δtf = 0.1
    u,v = unpack(sv,x)
    δu,δv = unpack(sv,δx)
    δun,δvn = DIVAnd_HFRadar.intertial_oscillation(Δtf,u,v,δu,δv)
    return pack(sv,(δun,δvn))[:,1]
end

function fun_adj(x,δx)
    Δtf = 0.1
    u,v = unpack(sv,x)
    δun,δvn = unpack(sv,δx)
    δu,δv = DIVAnd_HFRadar.intertial_oscillation_adj(Δtf,u,v,δun,δvn)
    return pack(sv,(δu,δv))[:,1]
end


#x1 = pack(sv,(randn((m-1,n)),randn((m,n-1))))[:,1]
#x2 = pack(sv,(randn((m-1,n)),randn((m,n-1))))[:,1]


x0 = randn(sv.n) # reference state, not used
x1 = randn(sv.n)
x2 = randn(sv.n)

#x1 = randn(2)
#x2 = randn(2)

#@show fun(x0,x2)
@test x1 ⋅ fun(x0,x2) ≈ fun_adj(x0,x1) ⋅ x2



