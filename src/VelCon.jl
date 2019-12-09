# velocity constrains

module VelCon

import DIVAnd
using Missings
using Distributed
using Statistics
using SparseArrays
using Dates
using LinearAlgebra

spdiagm(v::AbstractVector) = sparse(Diagonal(v))

include("stagger_mask.jl")
include("stagger_r2u.jl")
include("stagger_r2v.jl")

include("stagger_u2r.jl")
include("stagger_v2r.jl")

include("inertial_oscillations.jl")
include("inertial_oscillations_geo.jl")

include("DIVAnd_hfradar.jl")

end
