import DIVAnd_hfradar
using NCDatasets
using DataArrays
import MAT
import PhysOcean

include("DIVAnd_hfradar_load.jl")


# # compute and remove mean
# robs2 = copy(robs_all)
# robs2[flagcv_all] .= NaN
# mr = PhysOcean.nanmean(robs2,3);
# robs_all = robs_all .- mr;

selection = :cv

uri = Array{Float64}((length(lonr),length(latr),length(timerange)))
vri = Array{Float64}((length(lonr),length(latr),length(timerange)))
ηri = Array{Float64}((length(lonr),length(latr),length(timerange)))

# @show DIVAnd_hfradar.cverr(
#                xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#                lonr,latr,timerange,
#                mask2d,h,(3e3,3e3,0.),(0.,0.,3600.),1e-4,-1,-1,-1,0.,1.; selection=selection)
# # ncenter = 50
# # sum(valid) = 2907
# #   9.982437 seconds (3.60 M allocations: 227.285 MiB, 0.93% gc time)
# #   1.525246 seconds (689.19 k allocations: 90.263 MiB, 1.76% gc time)
# #   0.109364 seconds (225.96 k allocations: 66.231 MiB, 10.27% gc time)

# # ncenter = 52
# # sum(valid) = 2989
# #   0.068493 seconds (194.84 k allocations: 58.380 MiB, 11.21% gc time)
# #   0.111266 seconds (230.67 k allocations: 67.620 MiB, 8.13% gc time)
# #   0.110308 seconds (229.66 k allocations: 67.164 MiB, 8.20% gc time)

# # status = (4000.0, 4000.0, 0.0, 0.0, 0.0, 3600.0, 0.0001, -1, -1, -1, 0.0, 1.0, 0.056014775004292944)
# # 0.056014775004292944

# # 2D bc

# @show DIVAnd_hfradar.cverr(
#                              xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#                              lonr,latr,timerange,
#                              mask2d,h,(3e3,3e3,0.),(0.,0.,3600.),1e-4,1e-2,-1,-1,0.,1.; selection=selection)
# # ncenter = 50
# # sum(valid) = 2907
# #   0.208934 seconds (191.04 k allocations: 57.911 MiB, 68.54% gc time)
# #   0.078801 seconds (226.96 k allocations: 67.240 MiB, 13.44% gc time)
# #   0.077911 seconds (225.96 k allocations: 66.230 MiB, 12.75% gc time)
# # ncenter = 52
# # sum(valid) = 2989
# #   0.063869 seconds (194.84 k allocations: 58.380 MiB, 8.60% gc time)
# #   0.110357 seconds (230.67 k allocations: 67.620 MiB, 7.85% gc time)
# #   0.109260 seconds (229.66 k allocations: 67.164 MiB, 9.69% gc time)
# # status = (4000.0, 4000.0, 0.0, 0.0, 0.0, 3600.0, 0.0001, 0.01, -1, -1, 0.0, 1.0, 0.056367175900273694)
# # 0.056367175900273694


# # 2D div

# @show DIVAnd_hfradar.cverr(
#                              xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#                              lonr,latr,timerange,
#                              mask2d,h,(3e3,3e3,0.),(0.,0.,3600.),1e-4,-1,4.9e8,-1,0.,1.; selection=selection)
# # ncenter = 50
# # sum(valid) = 2907
# #   0.060722 seconds (191.04 k allocations: 57.911 MiB)
# #   0.102319 seconds (226.96 k allocations: 67.240 MiB)
# #   0.102117 seconds (225.96 k allocations: 66.230 MiB)
# # ncenter = 52
# # sum(valid) = 2989
# #   0.061728 seconds (194.84 k allocations: 58.380 MiB)
# #   0.106316 seconds (230.67 k allocations: 67.620 MiB)
# #   0.104222 seconds (229.66 k allocations: 67.164 MiB)
# # status = (4000.0, 4000.0, 0.0, 0.0, 0.0, 3600.0, 0.0001, -1, 4.9e8, -1, 0.0, 1.0, 0.05590168684002164)
# # 0.05590168684002164


# # 3D
# @show DIVAnd_hfradar.cverr(
#                       xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#                       lonr,latr,timerange,
#                       mask2d,h,(3e3,3e3,60*60 ),(0.,0.,3600.),1e-1,-1,-1,-1,0.,1.; selection=selection)
# # ncenter = 50
# # sum(valid) = 2907
# #   0.201958 seconds (191.04 k allocations: 57.911 MiB, 72.29% gc time)
# #   0.160785 seconds (286.75 k allocations: 113.146 MiB, 8.62% gc time)
# #   0.157136 seconds (285.28 k allocations: 111.792 MiB, 8.14% gc time)

# # ncenter = 52
# # sum(valid) = 2989
# #   0.067695 seconds (194.84 k allocations: 58.380 MiB, 12.87% gc time)
# #   0.160392 seconds (290.45 k allocations: 113.560 MiB, 9.12% gc time)
# #   0.158170 seconds (288.99 k allocations: 112.636 MiB, 8.05% gc time)


# # status = (4000.0, 4000.0, 3600, 0.0, 0.0, 3600.0, 0.1, -1, -1, -1, 0.0, 1.0, 0.0455705382844484)
# # 0.0455705382844484


# # Coriolis

# @show DIVAnd_hfradar.cverr(
#                xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#                lonr,latr,timerange,
#                mask2d,h,(3e3,3e3,0.),(0.,0.,24*10*3600.),1e-4,-1,-1,1e-5,0.,10.; selection=selection)
# # ncenter = 50
# # sum(valid) = 2907
# #   0.214497 seconds (190.70 k allocations: 58.906 MiB, 72.99% gc time)
# #   0.111810 seconds (226.96 k allocations: 67.240 MiB, 8.53% gc time)
# #   0.108661 seconds (225.96 k allocations: 66.230 MiB, 8.45% gc time)
# # success = true
# # ncenter = 52
# # sum(valid) = 2989
# #   0.070153 seconds (194.45 k allocations: 59.455 MiB, 13.95% gc time)
# #   0.112706 seconds (230.67 k allocations: 67.620 MiB, 11.70% gc time)
# #   0.111769 seconds (229.66 k allocations: 67.164 MiB, 10.47% gc time)
# # success = true
# # status = (4000.0, 4000.0, 0.0, 0.0, 0.0, 864000.0, 0.0001, -1, -1, 1.0e-5, 0.0, 10.0, 0.04384261562938635)
# # 0.04384261562938635


# Coriolis + geo


lenxy = 2889.25
    
#@show DIVAnd_hfradar.cverr(
#               xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#               lonr,latr,timerange,
#    mask2d,h,(3e3,3e3,0.),(0,0,24*10*3600.),1e-4,-1,-1,1e-5,9.81,7.; selection=selection)




# cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
#     xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
#     lonr,latr,timerange,
#     mask2d,htot,(lenxy,lenxy,0.),(0.,0.,24*10*3600),x[1],-1,-1,x[2],g_barotropic,x[3]; selection=selection,u = uri, v = vri, η = ηri)

# @show cverr3D_Coriolis_geo([5.981e-05, 4.347e-05, 1.138])
# 0.04835521434506809

cverr3D_Coriolis_geo2(x) = DIVAnd_hfradar.cverr(
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,htot,(lenxy,lenxy,0.),(0.,0.,24*5*3600),x[1],-1,-1,x[2],g_barotropic,x[3]; selection=selection)

@show cverr3D_Coriolis_geo2([5.981e-05, 4.347e-05, 1.138])


cverr3D_Coriolis_geo3(x) = DIVAnd_hfradar.cverr(
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,htot,(lenxy,lenxy,0.),(0.,0.,24*15*3600),x[1],-1,-1,x[2],g_barotropic,x[3]; selection=selection)

@show cverr3D_Coriolis_geo3([5.981e-05, 4.347e-05, 1.138])


# ncenter = 50
# sum(valid) = 2907
#   0.680822 seconds (300.53 k allocations: 64.129 MiB, 22.65% gc time)
#   0.110896 seconds (226.96 k allocations: 67.240 MiB, 8.73% gc time)
#   0.109274 seconds (225.96 k allocations: 66.230 MiB, 8.87% gc time)
# success = false
# ncenter = 52
# sum(valid) = 2989
#   0.066702 seconds (194.46 k allocations: 59.456 MiB, 8.41% gc time)
#   0.114203 seconds (230.67 k allocations: 67.620 MiB, 10.44% gc time)
#   0.113363 seconds (229.66 k allocations: 67.164 MiB, 8.17% gc time)
# success = false
# status = (4000.0, 4000.0, 0.0, 0, 0, 864000.0, 0.0001, -1, -1, 1.0e-5, 9.81, 7.0, 0.04188440835076011)
# 0.04188440835076011
