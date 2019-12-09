@everywhere using DIVAnd
@everywhere using GeoMapping
@everywhere import VelCon
@everywhere using NCDatasets
@everywhere using DataArrays
@everywhere using BlackBoxOptim
@everywhere import MAT
import PhysOcean


@everywhere include("DIVAnd_hfradar_load.jl")


# compute and remove mean
robs2 = copy(robs_all)
robs2[flagcv_all] .= NaN
mr = PhysOcean.nanmean(robs2,3);
robs_all = robs_all .- mr;

const postfix = get(ENV,"SLURM_JOBID","")
const cases = split(get(ENV,"SLURM_JOB_NAME",""),",")

#@show @__LINE__,@__FILE__

selection = :cv

#lenxy = 3000.
lenxy = 2889.25

if "2D_len" in cases
    # optmize lenxy and epsilon2
    cverr2D(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D;
                     SearchRange = [(2e3,20e3),
                                    (1e-5,1e-3)],
                     NumDimensions = 2, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr2D_len$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end

if "2D" in cases
    # optmize lenxy and epsilon2
    cverr2D(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(lenxy,lenxy,0.),(0.,0.,3600.),x[1],-1,-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D;
                     SearchRange = [(1e-5,1e-3)],
                     NumDimensions = 1, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr2D$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "2D_bc" in cases
    # optmize lenxy and epsilon2  bc
    cverr2D_bc(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(lenxy,lenxy,0.),(0.,0.,3600.),x[1],x[2],-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D_bc;
                     SearchRange = [(1e-5,1e-4),
                                    (1e-2,10.)],
                     NumDimensions = 2, TargetFitness=1e-4)

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
    cverr2D_div(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(lenxy,lenxy,0.),(0.,0.,3600.),x[1],-1,x[2],-1,0.,1.; selection=selection)

    res = bboptimize(cverr2D_div;
                     SearchRange = [(1e-5,1e-4),
                                    (1e8,1e9)],
                     NumDimensions = 2, TargetFitness=1e-4)

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
    cverr3D(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(lenxy,lenxy,x[1]),(0.,0.,3600.),x[2],-1,-1,-1,0.,1.; selection=selection)

    res = bboptimize(cverr3D;
                     SearchRange = [(0.5*60*60,2*60*60.),
                                    (1e-2,1.),
                                    ],
                     NumDimensions = 2, TargetFitness=1e-4)

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
    cverr3D_Coriolis(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(lenxy,lenxy,0.),(0.,0.,24*10*3600.),x[1],-1,-1,x[2],0.,1.; selection=selection)

    res = bboptimize(cverr3D_Coriolis;
                     SearchRange = [(1e-5,1e-3),
                                    (1e-6,1e-4),
                                    ],
                     NumDimensions = 2, TargetFitness=1e-4)

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
    cverr3D_Coriolis_geo(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(lenxy,lenxy,0.),(0.,0.,24*10*3600),x[1],-1,-1,x[2],9.81,x[3]; selection=selection)

    res = bboptimize(cverr3D_Coriolis_geo;
                     SearchRange = [(1e-5,1e-3),
                                    (1e-6,1e-4),
                                    (1.,50.),
                                    ],
                     NumDimensions = 3, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr3D_Coriolis_geo$(postfix).json","w") do f
        JSON.print(f,
                   Dict("xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end

@show "done"
#processall()
