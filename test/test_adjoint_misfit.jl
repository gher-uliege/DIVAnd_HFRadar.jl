using DIVAnd_hfradar
using DIVAnd
using Test
using LinearAlgebra


const m = 10
const n = 11
const nmax = 3
sv = statevector((trues((m-1,n,nmax)),trues((m,n-1,nmax))))

#x1 = pack(sv,(randn((m-1,n)),randn((m,n-1))))[:,1]
#x2 = pack(sv,(randn((m-1,n)),randn((m,n-1))))[:,1]

x0 = zeros(sv.n) # not used
δx1 = randn(sv.n)
δx2 = randn(sv.n)

#x1 = randn(2)
#x2 = randn(2)

t = 1:nmax;
f = 0.1;

@test δx1 ⋅ DIVAnd_hfradar.misfit_intertial_oscillation(sv,f,t,x0,δx2) ≈
    DIVAnd_hfradar.misfit_intertial_oscillation_adj(sv,f,t,x0,δx1) ⋅ δx2



const sv2 = statevector((trues((m,n,nmax)),trues((m-1,n,nmax)),trues((m,n-1,nmax))))

mask = trues(m,n)
h = ones(m,n)
pm = ones(m,n)
pn = ones(m,n)

t = 1:nmax;
f = 0.1;
g = 9.81;

grid = DIVAnd_hfradar.CGrid(mask,h,(pm,pn))
config = DIVAnd_hfradar.Config(grid,f,g)

#x1 = pack(sv,(randn((m-1,n)),randn((m,n-1))))[:,1]
#x2 = pack(sv,(randn((m-1,n)),randn((m,n-1))))[:,1]

x0 = zeros(sv2.n) # not used
δx1 = randn(sv2.n)
δx2 = randn(sv2.n)

#x1 = randn(2)
#x2 = randn(2)


@test δx1 ⋅ DIVAnd_hfradar.misfit_intertial_oscillation_geo(sv2,config,t,x0,δx2) ≈
    DIVAnd_hfradar.misfit_intertial_oscillation_geo_adj(sv2,config,t,x0,δx1) ⋅ δx2

#end
