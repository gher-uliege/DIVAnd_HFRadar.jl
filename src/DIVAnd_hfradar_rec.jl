import VelCon
using NCDatasets
import MAT
import PhysOcean
using GeoMapping
using JLD
using Glob
using Statistics
using DelimitedFiles

include("DIVAnd_hfradar_save.jl")
include("DIVAnd_hfradar_load.jl")
include("DIVAnd_hfradar_load_optim.jl")

# # compute and remove mean
# robs2 = copy(robs_all)
# robs2[flagcv_all] .= NaN
# mr = PhysOcean.nanmean(robs2,3);
# robs_all = robs_all .- mr;

selection = :all
#selection = :debug

postfix = "rg50"
postfix = "rg"
mkpath(expanduser("~/tmp/HFRadar-Ibiza/$(postfix)"))

uri = Array{Float64,3}(undef,(length(lonr),length(latr),length(timerange)))
vri = Array{Float64,3}(undef,(length(lonr),length(latr),length(timerange)))
ηri = Array{Float64,3}(undef,(length(lonr),length(latr),length(timerange)))


#cases = ["2D","2D_bc","2D_div","3D","3D_Coriolis","3D_Coriolis_geo","3D_Coriolis_geo_gp"]

cases = ["3D_Coriolis_geo"]
for c in cases

    RMS,param = load_optim_param(c,postfix)

    len = (param["lenxy"],param["lenxy"],param["lent"])
    lenη = (0.,0.,param["lenetat"])

    eps2 = param["eps2"]
    eps2_boundary_constrain = param["eps2_boundary_constrain"]
    eps2_Coriolis_constrain = param["eps2_Coriolis_constrain"]
    eps2_div_constrain = param["eps2_div_constrain"]
    g = param["g"]
    ratio = param["ratio"]


    @time @show VelCon.cverr(xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
                   lonr,latr,timerange,
                   mask2d,htot,
                   len,lenη,
                   eps2,
                   eps2_boundary_constrain,
                   eps2_div_constrain,
                   eps2_Coriolis_constrain,
                   g,
                   ratio;
                   selection=selection,u = uri, v = vri, η = ηri)


    fname = expanduser("~/tmp/HFRadar-Ibiza/$(postfix)/$(c).nc")

    DIVAnd_hfradar_save(fname,lonr,latr,timerange,uri,vri,ηri)
end
