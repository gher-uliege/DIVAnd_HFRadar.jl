import DIVAnd_hfradar
using NCDatasets
import MAT
import PhysOcean
using GeoMapping
#using JLD
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
#selection = :debug
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
cases = ["3D_Coriolis_geo_eta_optim2"]
cases = ["3D"]

c = cases[1]
#for c in cases
#=

    RMS,param = DIVAnd_hfradar.load_optim_param(c,postfix)

    @show RMS,param
    len = (param["lenxy"],param["lenxy"],param["lent"])
    lenη = (0.,0.,param["lenetat"])

    eps2 = param["eps2"]
    eps2_boundary_constraint = param["eps2_boundary_constraint"]
    eps2_Coriolis_constraint = param["eps2_Coriolis_constraint"]
    eps2_div_constraint = param["eps2_div_constraint"]
    g = param["g"]
    ratio = param["ratio"]

    @time @show DIVAnd_hfradar.cverr(xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
                   lonr,latr,timerange,
                   mask2d,htot,
                   len,lenη,
                   eps2,
                   eps2_boundary_constraint,
                   eps2_div_constraint,
                   eps2_Coriolis_constraint,
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


#=
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


=#


docv = true
docv = false

if !docv
    @info "disable cross-validation; use all data"
    flagcv_all .= false
end


#=
if "3D_Coriolis_geo_eta_optim2" in cases
    Δn = 1
    Δn = 2
    Δn = 3
    Δn = 4
    Δn = 5
    Δn = parse(Int,ENV["DELTA_N"])
    @show Δn

    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,(x[1],x[1],0.),(0*x[2],0*x[2],24*60*60*x[3]),x[4],-1,-1,x[5],g_barotropic,x[6]; selection=selection,
        u = uri, v = vri, η = ηri,
        Δn = Δn,
    )

end


    xopt = [16291.007298025252, 14603.129071943038, 7.906700502571824, 0.011196566797111344, 0.007994164838883752, 13.598806798689564]
    xopt = [9125.976646579633,                  0,  11.91390102142871,0.024069234801430514,0.012465245160690287,6.564149106945681]

    @show cverr3D_Coriolis_geo(xopt)


    #@show cverr3D_Coriolis_geo([5.981e-05 * 3600, 4.347e-05  * 3600, 1.138])


    fname = expanduser("~/tmp/HFRadar-Ibiza/$(postfix)/$(c)-Deltan$(Δn)-rerun-newopt.nc")

    DIVAnd_hfradar.DIVAnd_hfradar_save(fname,lonr,latr,timerange,uri,vri,ηri)
#end


=#


Δn = 1

    cverr3D(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],x[2]),(0.,0.,3600.),x[3],-1,-1,-1,0.,1.; selection=selection,
        u = uri, v = vri, η = ηri,
        Δn = Δn,
    )
xopt = [2134.307509439117,3526.195062616695,0.15101748832935727]
    @show cverr3D(xopt)

fname = expanduser("~/tmp/HFRadar-Ibiza/$(postfix)/$(c)-Deltan$(Δn)-rerun.nc")

    DIVAnd_hfradar.DIVAnd_hfradar_save(fname,lonr,latr,timerange,uri,vri,ηri)
