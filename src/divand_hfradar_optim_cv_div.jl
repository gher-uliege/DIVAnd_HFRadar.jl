@everywhere using DIVAnd
@everywhere using GeoMapping
@everywhere import VelCon
@everywhere using NCDatasets
@everywhere using DataArrays
@everywhere using BlackBoxOptim
@everywhere import MAT

# bearing β: angle at radar station (*) between North a measuring point (+) counted clockwise
#
# direction α: angle at measuring point between North and vector pointing to the radar station counted clockwise
#
#              ↑ /
#              |/ 
#          ↑   +--→ current vector (u,v)
#   North  |  / measurent point
#          |β/
#          |/
#          *
#    radar station
#           
# direction α = bearing β + 180
#
# u zonal component, v meridional component
# u = r * sin(α)
# v = r * cos(α)  
#
# r = u * sin(α) + v * cos(α) 
#
# 
# u = -r * sin(β)
# v = -r * cos(β)
#
# r, u, v, direction and β consistent with the CODAR convention of the ruv files

#%%   Longitude   Latitude    U comp   V comp  VectorFlag    Spatial    Temporal     Velocity    Velocity  Spatial  Temporal X Distance  Y Distance   Range   Bearing   Velocity  Direction   Spectra
#%%     (deg)       (deg)     (cm/s)   (cm/s)  (GridCode)    Quality     Quality     Maximum     Minimum    Count    Count      (km)        (km)       (km)    (True)    (cm/s)     (True)    RngCell
#     1.3721884  38.6737279   -5.032    2.904          0       1.085       2.478      -5.267      -6.352       2        4      -1.4412      0.8321    1.6642   300.0     -5.810     120.0         1
#     1.3409888  38.6013079    4.105    7.118          0       2.241       1.506      17.520       6.669       5        7      -4.1605     -7.2062    8.3210   210.0      8.217      30.0         5


# Example
#  u = -5.032 = -5.810 * sin(120.0 * pi/180)
#  v =  2.904 = -5.810 * cos(120.0 * pi/180)
#
#  u = 4.105 = 8.217 * sin(30.0 * pi/180)
#  v = 7.118 = 8.217 * cos(30.0 * pi/180)

# r is positive if velocity is pointing *towards* the site.
# e.g. for example u = 1, v = 0 (water going from the west to the east) and bearing = 90°, direction = 90° + 180° = 270°
#
#       N ↑
#         |
#         *-----→ current vector (u,v) = (1,0)
#    radar station
#
# r = 1 * sin(270°) = -1. The radial velocity r is negative, thus the currents are away from the site.

# CODAR manual
# "the sign of the velocity indicates whether the velocity is towards (positive) or away (negative) from the site" (page 12)
# http://support.codar.com/Technicians_Information_Page_for_SeaSondes/Docs/GuidesToFileFormats/File_LonLatUV_RDL_TOT_ELP.pdf
# https://web.archive.org/web/20170801200228/http://support.codar.com/Technicians_Information_Page_for_SeaSondes/Docs/GuidesToFileFormats/File_LonLatUV_RDL_TOT_ELP.pdf




#converttime(dt) = Int(dt-timeorigin)/1000

@everywhere include("DIVAnd_hfradar_load.jl")




if false
    for n = 1:length(ntimes)
        figure()

        dt = (Int(timerange[n] - timeorigin))/1000
        hfradar_plot(xi[:,:,n],yi[:,:,n],dt,uri[:,:,n],vri[:,:,n],
                     xobs,yobs,tobs,robs,directionobs,sitenames_,siteorigin)
    end

end


postfix = get(ENV,"SLURM_JOBID","")

#@show @__LINE__,@__FILE__

             
cverr2D(x) = VelCon.cverr(
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,-1,-1,0.,1.)


@show cverr2D([3e3,1e-4])


if true
    # optmize lenxy and epsilon2 div
    cverr2D_div(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,x[3],-1,0.,1.)

    res = bboptimize(cverr2D_div;
                     SearchRange = [(1e3, 6e3),
                                    (1e-5,1e-3),
                                    (1e-5,100.)],
                     NumDimensions = 3, TargetFitness=1e-4)
    
    xopt = best_candidate(res)
    @show xopt
    
    open("cverr2D_div$(postfix).json","w") do f
        JSON.print(f,
                   Dict("lenxy" => xopt[1],
                        "epsilon2" => xopt[2],
                        "eps2_div_constrain" => xopt[3],
                        "xopt" => xopt,
                        "best_fitness" => best_fitness(res)));    
    end
end

#processall()
