@everywhere using DIVAnd
@everywhere using GeoMapping
@everywhere import VelCon
@everywhere using NCDatasets
@everywhere using DataArrays
@everywhere using BlackBoxOptim
@everywhere import MAT


@everywhere include("DIVAnd_hfradar_load.jl")


const postfix = get(ENV,"SLURM_JOBID","")
const cases = split(get(ENV,"SLURM_JOB_NAME",""),",")

#@show @__LINE__,@__FILE__

             

if "2D" in cases
    # optmize lenxy and epsilon2
    cverr2D(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,-1,-1,0.,1.)

    res = bboptimize(cverr2D; 
#                     SearchRange = [(2e3, 6e3),
#                                    (1e-5,1e-3)],
                     SearchRange = [(4e3, 8e3),
                                    (1e-5,1e-3)],
                     NumDimensions = 2, TargetFitness=1e-4)

    xopt = best_candidate(res)
    @show xopt

    open("cverr2D$(postfix).json","w") do f
        JSON.print(f,
                   Dict("lenxy" => xopt[1],
                        "epsilon2" => xopt[2],
                        "xopt" => xopt,
                        "best_fitness" => best_fitness(res)));    
    end
end


if "2D_bc" in cases
    # optmize lenxy and epsilon2  bc
    cverr2D_bc(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],x[3],-1,-1,0.,1.)

    res = bboptimize(cverr2D_bc;
                     SearchRange = [(4e3, 8e3),
                                    (1e-5,1e-4),
                                    (1e-5,100.)],
                     NumDimensions = 3, TargetFitness=1e-4)
    
    xopt = best_candidate(res)
    @show xopt
    
    open("cverr2D_bc$(postfix).json","w") do f
        JSON.print(f,
                   Dict("lenxy" => xopt[1],
                        "epsilon2" => xopt[2],
                        "eps2_boundary_constrain" => xopt[3],
                        "xopt" => xopt,
                        "best_fitness" => best_fitness(res)));
    end
end


if "2D_div" in cases
    # optmize lenxy and epsilon2 div
    cverr2D_div(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,x[3],-1,0.,1.)

    res = bboptimize(cverr2D_div;
                     SearchRange = [(3e3, 8e3),
                                    (1e-5,1e-4),
                                    (1e-5,1e4)],
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


if "3D" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],x[2]),(0.,0.,3600.),x[3],-1,-1,-1,0.,1.)

    res = bboptimize(cverr3D;
                     SearchRange = [(3e3, 8e3),
                                    (0.,2*60*60.),
                                    (1e-5,1e-4),
                                    ],
                     NumDimensions = 3, TargetFitness=1e-4)
    
    xopt = best_candidate(res)
    @show xopt
    
    open("cverr3D$(postfix).json","w") do f
        JSON.print(f,
                   Dict("lenxy" => xopt[1],
                        "lent" => xopt[2],
                        "epsilon2" => xopt[3],
                        "xopt" => xopt,
                        "best_fitness" => best_fitness(res)));    
    end
end


if "3D_Coriolis" in cases
    # optmize lenxy and epsilon2 div and lent
    cverr3D_Coriolis(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,3600.),x[2],-1,-1,x[3],0.,1.)

    res = bboptimize(cverr3D_Coriolis;
                     SearchRange = [(3e3, 8e3),
                                    (1e-5,1e-4),
                                    (1e8,1e10),
                                    ],
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
    cverr3D_Coriolis_geo(x) = VelCon.cverr(
        xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
        lonr,latr,timerange,
        mask2d,h,(x[1],x[1],0.),(0.,0.,x[2]),x[3],-1,-1,x[4],9.81,x[5])

    res = bboptimize(cverr3D_Coriolis_geo;
                     SearchRange = [(3e3, 8e3),
                                    (1.*60.,3*60*60.),
                                    (1e-5,1e-4),
                                    (1e8,1e10),
                                    (1.,200.),
                                    ],
                     NumDimensions = 5, TargetFitness=1e-4)
    
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
