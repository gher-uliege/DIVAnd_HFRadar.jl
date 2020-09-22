using Distributed
@everywhere using DIVAnd
@everywhere using GeoMapping
@everywhere import DIVAnd_hfradar
@everywhere using NCDatasets
@everywhere using Statistics
@everywhere using BlackBoxOptim
@everywhere import MAT
import PhysOcean
import JSON

@everywhere include(joinpath(dirname(@__FILE__),"DIVAnd_hfradar_load.jl"))
#include("DIVAnd_hfradar_load.jl")


const postfix = get(ENV,"SLURM_JOBID","")
const cases = split(get(ENV,"SLURM_JOB_NAME",""),",")

#@show @__LINE__,@__FILE__

selection = :cv

#lenxy = 3000.
#lenxy = 2889.25

if "2D" in cases
    # optmize lenxy and epsilon2
    cverr2D(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,
        (x[1],x[1],0.),
        (0.,0.,3600.),
        x[2],-1,-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D;
                     SearchRange = [
                         (2e3,50e3),
                         (1e-5,1e-1),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 2, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr2D_len$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "2D_bc" in cases
    # optmize lenxy and epsilon2  bc
    cverr2D_bc(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,
        (x[1],x[1],0.),
        (0.,0.,3600.),
        x[2],x[3],-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D_bc;
                     SearchRange = [
                         (2e3,50e3),
                         (1e-5,1e-1),
                         (1e-2,10.),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 3, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr2D_bc$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "2D_div" in cases
    # optmize lenxy and epsilon2 div
    cverr2D_div(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,x[3],-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D_div;
                     SearchRange = [
                         (2e3,80e3),
                         (1e-5,1e-1),
                         (1e9,1e12),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 3, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr2D_div$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "3D" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],x[2]),(0.,0.,3600.),x[3],-1,-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr3D;
                     SearchRange = [
                         (2e3,50e3),
                         (0.5*60*60,60*60.),
                         (1e-2,1.),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 3, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end



if "3D_Coriolis" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,24*10*3600.),x[2],-1,-1,x[3],0.,1.; selection=selection)

    res = bboptimize(cverr3D_Coriolis;
                     SearchRange = [
                         (2e3,50e3),
                         (1e-5,1e-1),
                         (1e-5,1e-2),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 3, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "3D_Coriolis_geo" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,(x[1],x[1],0.),(0.,0.,24*10*3600.),x[2],-1,-1,x[3],g_barotropic,x[4]; selection=selection)

    res = bboptimize(cverr3D_Coriolis_geo;
                     SearchRange = [
                         (2e3,50e3),
                         (1e-3,1e-1),
                         #(1e-6,1e-4),
                         (1e-3,1e-1),
                         (1.,50.),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 4, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis_geo$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end

if "3D_Coriolis_geo_eta_optim2" in cases

    Δn = 1
    Δn = 2
    Δn = 3
    @show Δn
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,
        #(x[1],x[1],0.),(0*x[2],0*x[2],24*60*60*x[3]),x[4],-1,-1,x[5],g_barotropic,x[6]; 
        (x[1],x[1],0.),(0.,0.,24*60*60*x[2]),x[3],-1,-1,x[4],g_barotropic,x[5]; 
        selection=selection, Δn = Δn)

    res = bboptimize(cverr3D_Coriolis_geo;
                     SearchRange = [
                         (2e3,50e3), # lenxy u,v
                         #(2e3,50e3),  # lenxy eta
                         (3.,12.),  # lent eta
                         (1e-3,1e-1), # epsilon obs
                         (1e-3,1e-1), # epsilon Coriolis_geo
                         (1.,50.), # ratio
                     ],
                     MaxSteps = 500,
                     #NumDimensions = 6,
                     NumDimensions = 5,
                     TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis_geo_Δn$(Δn)_$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end



if "3D_Coriolis_geo_eta_optim_lenxy" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,(x[1],x[1],0.),(x[2],x[2],24*60*60*x[3]),x[4],-1,-1,x[5],g_barotropic,x[6]; selection=selection)

    res = bboptimize(cverr3D_Coriolis_geo;
                     SearchRange = [
                         (2e3,50e3), # lenxy u,v
                         (2e3,50e3),  # lenxy eta
                         (3.,12.),  # lent eta
                         (1e-3,1e-1), # epsilon obs
                         (1e-3,1e-1), # epsilon Coriolis_geo
                         (1.,50.), # ratio
                     ],
                     MaxSteps = 500,
                     NumDimensions = 6, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis_geo$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "3D_Coriolis_geo_eta_optim_lenxy_lent" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis_geo(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,htot,(x[1],x[1],x[2]),(x[3],x[3],24*60*60*x[4]),x[6],-1,-1,x[6],g_barotropic,x[7]; selection=selection)

    res = bboptimize(cverr3D_Coriolis_geo;
                     SearchRange = [
                         (2e3,50e3), # lenxy u,v
                         (0.5*60*60,60*60.), # lent
                         (2e3,50e3),  # lenxy eta
                         (3.,12.),  # lent eta
                         (1e-3,1e-1), # epsilon obs
                         (1e-3,1e-1), # epsilon Coriolis_geo
                         (1.,50.), # ratio
                     ],
                     MaxSteps = 500,
                     NumDimensions = 7, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis_geo$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end





if "3D_Coriolis_geo_gp" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis_geo_gp(x) = DIVAnd_hfradar.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,24*10*3600.),x[2],-1,-1,x[3],g_baroclinic,x[4]; selection=selection)

    res = bboptimize(cverr3D_Coriolis_geo_gp;
                     SearchRange = [
                         (2e3,50e3),
                         (1e-5,1e-1),
                         #(1e-5,1e-2),
                         (1e-3,1e-1),
                         (1e-2,2.),
                     ],
                     MaxSteps = 500,
                     NumDimensions = 4, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis_geo_gp$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end



print("done")
#processall()
