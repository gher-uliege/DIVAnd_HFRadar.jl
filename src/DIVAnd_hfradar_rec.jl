import DIVAnd_hfradar
using NCDatasets
import MAT
import PhysOcean
using GeoMapping
using JLD
using Glob
using Statistics
using DelimitedFiles

include("DIVAnd_hfradar_load.jl")

# # compute and remove mean
# robs2 = copy(robs_all)
# robs2[flagcv_all] .= NaN
# mr = PhysOcean.nanmean(robs2,3);
# robs_all = robs_all .- mr;

selection = :all
selection = :debug
#selection = :cv

postfix = "rg50"
#postfix = "rg"
mkpath(expanduser("~/tmp/HFRadar-Ibiza/$(postfix)"))

lenxy = 2889.25
eps2 = 0.0001161

uri = Array{Float64,3}(undef,(length(lonr),length(latr),length(timerange)))
vri = Array{Float64,3}(undef,(length(lonr),length(latr),length(timerange)))
ηri = Array{Float64,3}(undef,(length(lonr),length(latr),length(timerange)))


#cases = ["2D","2D_bc","2D_div","3D","3D_Coriolis","3D_Coriolis_geo","3D_Coriolis_geo_gp"]

cases = ["3D_Coriolis_geo"]
c = cases[1]
#for c in cases
#=

    RMS,param = DIVAnd_hfradar.load_optim_param(c,postfix)

    @show RMS,param
    len = (param["lenxy"],param["lenxy"],param["lent"])
    lenη = (0.,0.,param["lenetat"])

    eps2 = param["eps2"]
    eps2_boundary_constrain = param["eps2_boundary_constrain"]
    eps2_Coriolis_constrain = param["eps2_Coriolis_constrain"]
    eps2_div_constrain = param["eps2_div_constrain"]
    g = param["g"]
    ratio = param["ratio"]

    @time @show DIVAnd_hfradar.cverr(xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
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
=#
#=
cverr2D(x) = DIVAnd_hfradar.cverr(
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,h,
    (lenxy,lenxy,0.),
    (0.,0.,3600.),
    x[1],-1,-1,-1,0.,1.;
    selection=selection,u = uri, v = vri, η = ηri)

@show cverr2D([eps2])
=#

#=
cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,
        (lenxy,lenxy,0.),(0.,0.,24*10*3600.),
        x[1],
        -1,
        -1,
        x[2],
        g_barotropic,
        x[3];
        selection=selection,u = uri, v = vri, η = ηri)

    @show cverr3D_Coriolis_geo([5.981e-05, 4.347e-05, 1.138])

=#

DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,
        (lenxy,lenxy,0.),(0.,0.,24*10*3600.),
#        5.981e-05 #= * 3600 =#,
        5.981e-05 * 3600,
        -1,
        -1,
#        4.347e-05 #=* 3600 =#,
        4.347e-05 * 3600,
        g_barotropic,
        1.138;
        selection=selection,u = uri, v = vri, η = ηri)





    #@show cverr3D_Coriolis_geo([5.981e-05 * 3600, 4.347e-05  * 3600, 1.138])


    fname = expanduser("~/tmp/HFRadar-Ibiza/$(postfix)/$(c)-rerun.nc")

    DIVAnd_hfradar.DIVAnd_hfradar_save(fname,lonr,latr,timerange,uri,vri,ηri)
#end
