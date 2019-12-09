import VelCon
using NCDatasets
using DataArrays
import MAT
import PhysOcean
using PyPlot
using OceanPlot
using GeoMapping
using JLD

include("DIVAnd_hfradar_save.jl")
include("DIVAnd_hfradar_load.jl")

# # compute and remove mean
# robs2 = copy(robs_all)
# robs2[flagcv_all] .= NaN
# mr = PhysOcean.nanmean(robs2,3);
# robs_all = robs_all .- mr;

selection = :cv

figdir = expanduser("~/Doc/DIVAnd_hfradar/Fig/")

selection = :debug
#selection = :all

lenxy = 2889.25
eps2 = 0.0001161


# 2D
uri = Array{Float64}((length(lonr),length(latr),length(timerange)))
vri = Array{Float64}((length(lonr),length(latr),length(timerange)))
ηri = Array{Float64}((length(lonr),length(latr),length(timerange)))

cverr2D(x) = VelCon.cverr(
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,h,(lenxy,lenxy,0.),(0.,0.,3600.),x[1],-1,-1,-1,0.,1.; selection=selection,u = uri, v = vri, η = ηri)

@show cverr2D([eps2])

fname = expanduser("~/tmp/HFRadar-Ibiza/2D.nc")
DIVAnd_hfradar_save(fname,lonr,latr,timerange,uri,vri,ηri)

include("hfradar_plot.jl");

xi,yi = DIVAnd.ndgrid(lonr,latr)
n = 52;
hfradar_plot(xi,yi,timerange[n],uri[:,:,n],vri[:,:,n],
             xobs_all,yobs_all,robs_all[:,:,n,:],directionobs_all[:,:,n,:],sitenames,siteorigin)

savefig("$(figdir)/DIVAnd_hfradar_2D.png")
savefig("$(figdir)/DIVAnd_hfradar_2D.svg")

# # 3D_Coriolis_geo

# cverr3D_Coriolis_geo(x) = VelCon.cverr(
#     xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#     lonr,latr,timerange,
#     mask2d,htot,(lenxy,lenxy,0.),(0.,0.,24*10*3600),x[1],-1,-1,x[2],g_barotropic,x[3]; selection=selection,u = uri, v = vri, η = ηri)

# @show cverr3D_Coriolis_geo([5.981e-05, 4.347e-05, 1.138])

# hfradar_plot(xi,yi,timerange[n],uri[:,:,n],vri[:,:,n],
#              xobs_all,yobs_all,robs_all[:,:,n,:],directionobs_all[:,:,n,:],sitenames,siteorigin)

# savefig("$(figdir)/DIVAnd_hfradar_3D_Coriolis_geo.png")
# savefig("$(figdir)/DIVAnd_hfradar_3D_Coriolis_geo.svg")
