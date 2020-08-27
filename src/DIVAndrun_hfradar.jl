using JLD

function spobsoper_radvel(sv,modelgrid,Xobs,varindexu,varindexv,direction)

    #@show size(sv.mask[varindexu])
    #@show size(modelgrid[varindexu])
    #@show size(Xobs[1])

    Hu,outu,outboxu = DIVAnd.sparse_interp_g(modelgrid[varindexu],sv.mask[varindexu],Xobs)
    Hv,outv,outboxv = DIVAnd.sparse_interp_g(modelgrid[varindexv],sv.mask[varindexv],Xobs)

    d = direction[:]*pi/180

    Hgridu = spzeros(length(Xobs[1]),sv.n)
    Hgridu[:,sv.ind[varindexu]+1:sv.ind[varindexu+1]] = spdiagm(sin.(d)) * Hu * DIVAnd.sparse_pack(sv.mask[varindexu])'

    Hgridv = spzeros(length(Xobs[1]),sv.n)
    Hgridv[:,sv.ind[varindexv]+1:sv.ind[varindexv+1]] = spdiagm(cos.(d)) * Hv * DIVAnd.sparse_pack(sv.mask[varindexv])'

    valid = .!outu .& .!outv
    H = DIVAnd.sparse_pack(valid) * (Hgridu + Hgridv)
    return H,valid
end

"""
Create a sparse matrix which extract all elements of a state vector correspond to a true value in masks.
masks is a tulple of boolean mask.
"""
function DIVAnd.sparse_pack(sv,masks)
    j = vcat([findall(@view masks[i][:]) .+ sv.ind[i]  for i in 1:length(masks)]...)
    return sparse(1:length(j),j,ones(size(j)),length(j),sv.n)
end


function cgrid(mask,h,xyi,pmn)
    # staggering

    # mask_u/v keeps also the velocity normal to the boundary so that we can impose its value to zero (or close to zero)

    mask_u, mask_v, mask_psi = stagger_mask(mask,|)

    h_u = stagger_r2u(h)
    h_v = stagger_r2v(h)

    xyi_u = ([stagger_r2u(xi) for xi in xyi]...,)
    pmn_u = ([stagger_r2u(pm) for pm in pmn]...,)

    xyi_v = ([stagger_r2v(xi) for xi in xyi]...,)
    pmn_v = ([stagger_r2v(pm) for pm in pmn]...,)


    return mask_u,mask_v,h_u,h_v,xyi_u,xyi_v,pmn_u,pmn_v
end


"""

    DIVAnd_hfradar(mask,h,pmn,xyi,xyobs,robs,directionobs,len,epsilon2;...)

HF Radar current analysis with DIVAnd and velocity contraints.

mask: true for sea points (false for land points) (3D-array)
h: depth in meters (3D-array)
pmn: inverse of the local resolution (tuple of three 3D-arrays)
xyi: coordinates of the analysis grid (tuple of three 3D-arrays)
xyobs: coordinates of the observations (tuple of three vectors)
robs: radial velocity (vector)
directionobs: angle α of the measured direction in degrees (vector) such that

```math
u\\_{obs} * sin(α) + v\\_{obs} * cos(α) ≈ robs
```

len,epsilon2
f in  s⁻¹


## Convention for the direction

bearing β: angle at radar station (*) between North a measuring point (+) counted clockwise
direction α: angle at measuring point between North and vector pointing to the radar station counted clockwise



                ↑ /
                |/
         ↑      +--→ current vector (u,v)
  North  |     / measurent point
         |    /
         |   /
         |  /
         |β/
         |/
         *
   radar station

direction α = bearing β + 180

u zonal component, v meridional component
u = r * sin(α)
v = r * cos(α)

r = u * sin(α) + v * cos(α)


u = -r * sin(β)
v = -r * cos(β)

For HF radar, r is positive if velocity is pointing *towards* the radar site.
r, u, v, direction and β consistent with the CODAR convention of the ruv files

"""
function DIVAndrun_hfradar(mask,h,pmn,xyi,xyobs,robs,directionobs,len,epsilon2;
                        eps2_boundary_constrain = -1,
                        eps2_div_constrain = -1,
                        eps2_Coriolis_constrain = -1,
                        f = 0.001,
                        residual = zeros(size(robs)),
                        g = 0., # no pressure (testing)
                        ratio = 100,
                        lenη = (000.0, 000.0, 24 * 60 * 60. * 10)
                        )



    sz = size(mask)

    if any([size(pm) != sz for pm in pmn])
        error("size of mask ($(DIVAnd.formatsize(sz))) does not match size of pmn ($(join([DIVAnd.formatsize(size(pm)) for pm in pmn],", ")))")
    end

    if any([size(xi) != sz for xi in xyi])
        error("size of mask ($(DIVAnd.formatsize(sz))) does not match size of xyi ($(join([DIVAnd.formatsize(size(xi)) for xi in xyi],", ")))")
    end

    if size(h) != sz
        error("size of mask ($(DIVAnd.formatsize(sz))) does not match size of h ($(DIVAnd.formatsize(size(h))))")
    end

    if any([size(xobs) != size(robs) for xobs in xyobs])
        error("size of robs ($(DIVAnd.formatsize(size(robs)))) does not match size of xyobs ($(join([DIVAnd.formatsize(size(xobs)) for xobs in xyobs],", ")))")
    end

    if size(robs) != size(directionobs)
        error("size of robs ($(DIVAnd.formatsize(size(robs)))) does not match size of directionobs ($(DIVAnd.formatsize(size(directionobs))))")
    end

    if isa(epsilon2,Vector)
        if size(robs) != size(epsilon2)
            error("size of robs ($(DIVAnd.formatsize(size(robs)))) does not match size of epsilon2 ($(DIVAnd.formatsize(size(epsilon2))))")
        end
    end

    # 3D
    grid = CGrid(mask[:,:,1],h[:,:,1],(pmn[1][:,:,1],pmn[2][:,:,1]))
    config = Config(grid,f,g)

    mask_u,mask_v,h_u,h_v,xyi_u,xyi_v,pmn_u,pmn_v = cgrid(mask,h,xyi,pmn)

    # boundary_psi is unused
    boundary_u, boundary_v, boundary_psi = stagger_mask(mask,xor)


    varindexu = 1
    varindexv = 2

    modelgrid = (xyi,xyi_u,xyi_v)

    sv = DIVAnd.statevector((mask,mask_u,mask_v));
    #@show sv.n
    #@show sum.(xyobs)
    H,valid = spobsoper_radvel(sv,modelgrid,xyobs,2,3,directionobs)
    #JLD.@save "/tmp/new_data.jld" H valid xyobs directionobs sv

    #o = JLD.load("/tmp/old_data.jld")
    #valid = o["valid"]
    #H = sparse(o["Hi"],o["Hj"],o["Hval"],sum(valid),sv.n)

    #@show sum(valid)
    #@show size(H)
    #@show length.(xyobs)
    #@show valid

    yo = robs[valid]

    #@show sum(mask), sum.(pmn), lenη

    @time fi_η,s_η = DIVAnd.DIVAndrun(mask,pmn,xyi,xyobs,robs,lenη,1., alphabc = 1);
    @time fi,s_u = DIVAnd.DIVAndrun(mask_u,pmn_u,xyi_u,xyobs,robs,len,1., alphabc = 1);
    @time fi,s_v = DIVAnd.DIVAndrun(mask_v,pmn_v,xyi_v,xyobs,robs,len,1., alphabc = 1);

    #@show sum(abs.(s_u.iB))
    iB = blockdiag(ratio * s_η.iB,s_u.iB,s_v.iB)

    #@show sum(iB)
    #@show sum(abs.(iB))
    #@show sum(iB.^2)
    # ignore cross-validation points in analysis

    if isa(epsilon2,Number)
        R = Diagonal(epsilon2 * ones(size(yo)))
    else
        R = Diagonal(epsilon2[valid])
    end

    # @show size(valid)
    # @show size(yo)
    # @show size(R)
    # @show size(H)

    iP = iB + H' * (R \ H)

    Pxa = H' * (R \ yo)
    #@show sum(Pxa.^2)

    if eps2_boundary_constrain != -1
        Hboundary = DIVAnd.sparse_pack(sv,(falses(size(mask)),boundary_u,boundary_v))
        Rboundary = eps2_boundary_constrain * Diagonal(ones(size(Hboundary,1)))
        yoboundary = zeros(size(Hboundary,1))

        iP += Hboundary' * (Rboundary \ Hboundary)
        Pxa +=  Hboundary' * (Rboundary \ yoboundary)
    end

    if eps2_div_constrain != -1
        # divergence
        #   ∂hu/∂x + ∂hv/∂y  ≈ 0
        # cost function
        #   (∂hu/∂x + ∂hv/∂y)² / ϵ²_div

        Iu = DIVAnd.sparse_diag(h_u[:] ./ pmn_u[2][:]);
        Iv = DIVAnd.sparse_diag(h_v[:] ./ pmn_v[1][:]);

        Nt = if length(sz) == 2
            1
        else
            sz[3]
        end

        TUy = DIVAnd.sparse_trim((sz[1]-1,sz[2],Nt),2);
        TVx = DIVAnd.sparse_trim((sz[1],sz[2]-1,Nt),1);
        DUx = DIVAnd.sparse_diff((sz[1]-1,sz[2]-2,Nt),1);
        DVy = DIVAnd.sparse_diff((sz[1]-2,sz[2]-1,Nt),2);

        Pu = DIVAnd.sparse_pack(mask_u);
        Pv = DIVAnd.sparse_pack(mask_v);

        # the first block is the surface elevation and it is not used
        # in the computation of the divergence
        UP = blockdiag(spzeros(0,sum(mask)),copy(Pu'),copy(Pv'));
        Hdiv = [DUx*TUy*Iu   DVy*TVx*Iv] * UP;
        Rdiv = eps2_div_constrain * Diagonal(ones(size(Hdiv,1)))

        #@show size(Hdiv)

        iP = iP + Hdiv' * (Rdiv \ Hdiv)
    end

    #@show @__LINE__,@__FILE__

    if eps2_Coriolis_constrain != -1
        #if false
        x0 = zeros(sv.n) # not used for a linear constrain
        ti = xyi[3]
        t = ti[1,1,:]

        function iPfun(x,iPx)
            iPx[:] = iP*x + (1/eps2_Coriolis_constrain) *
                misfit_intertial_oscillation_geo_adj(
                    sv,config,t,x0,misfit_intertial_oscillation_geo(
                        sv,config,t,x0,x))
        end

        fi,success,niter = DIVAnd.conjugategradient(
            iPfun,Pxa,tol=1e-6,maxit = 50000,
            #        progress = DIVAnd.cgprogress
        )

        @show success
    else
        fi = iP \ Pxa
    end

    ηi,ui,vi = DIVAnd.statevector_unpack(sv,fi,NaN)

    uri = stagger_u2r(ui)
    vri = stagger_v2r(vi)

    residual[:] .= NaN
    residual[valid] = H*fi - yo

    #@show sum(uri[:,:,2].^2)
    #display(uri[:,:,2])

    return uri,vri,ηi
end




function cv_rec(
    ncenter,Δn,
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,h,
    len,lenη,eps2,
    eps2_boundary_constrain,
    eps2_div_constrain,
    eps2_Coriolis_constrain,
    g,
    ratio
)

    timeorigin = DateTime(2000,1,1)
    converttime(dt) = Float64(Dates.Second(dt-timeorigin).value)

    Ω = 7.2921e-5 # rad/s
    f = 2*Ω*sin(mean(latr) * π /180)

    imax = size(xobs_all,1)
    jmax = size(xobs_all,2)
    nsites = length(sitenames)
    #@show @__LINE__,@__FILE__

    #ncenter = ncv[l]
    @show ncenter
    ntimes = max(ncenter-Δn,1):min(ncenter+Δn,length(timerange))

    # select valid points within time range defined by ntimes
    valid = .!isnan.(robs_all[:,:,ntimes,:])

    robs = robs_all[:,:,ntimes,:][valid]
    directionobs = directionobs_all[:,:,ntimes,:][valid]

    forcv = flagcv_all[:,:,ntimes,:][valid]

    xobs = repeat(
        reshape(xobs_all,(imax,jmax,1,size(xobs_all,3)));
        inner = (1,1,length(ntimes),1))[valid]

    yobs = repeat(
        reshape(yobs_all,(imax,jmax,1,size(xobs_all,3)));
        inner = (1,1,length(ntimes),1))[valid]

    #sitenames_ = sitenames[ind2sub(size(valid),findall(valid))[4]]


    tobs = repeat(reshape(converttime.(timerange[ntimes]),(1,1,length(ntimes),1)),inner = (imax,jmax,1,nsites))[valid];

    residual = fill(NaN,size(robs))

    if length(timerange) == 1
        xi,yi = DIVAnd.ndgrid(lonr,latr)

        pm = ones(size(xi)) / (DIVAnd.deg2m(1) * (xi[2,1]-xi[1,1]));
        pn = ones(size(xi)) / (DIVAnd.deg2m(1) * cos(mean(yi) * pi/180) * (yi[1,2]-yi[1,1]));

        #len = ntuple(i -> lenxy,length(sz))

        pmn = (pm,pn)
        xyi = (xi,yi)
        xyobs = (xobs,yobs)

        epsilon2 = fill(eps2,size(robs))
    else
        mask = repeat(mask2d,inner=(1,1,length(ntimes)))
        h3d = repeat(h,inner=(1,1,length(ntimes)))
        xi,yi,ti = DIVAnd.ndgrid(lonr,latr,converttime.(timerange[ntimes]))

        pm = ones(size(xi)) / (DIVAnd.deg2m(1) * (xi[2,1,1]-xi[1,1,1]));
        pn = ones(size(xi)) / (DIVAnd.deg2m(1) * cos(mean(yi) * pi/180) * (yi[1,2,1]-yi[1,1,1]));
        po = ones(size(xi)) / (ti[1,1,2]-ti[1,1,1]);

        # correlation length in meters
        #len = (lenxy,lenxy,lent)

        pmn = (pm,pn,po)
        xyi = (xi,yi,ti)
        xyobs = (xobs,yobs,tobs)

        epsilon2 = fill(eps2,size(robs))
    end


    sz = size(mask)
    #h = ones(sz)


    epsilon2[forcv] .= Inf


    uri_temp,vri_temp,ηri_temp = DIVAndrun_hfradar(
        mask,h3d,pmn,xyi,xyobs,robs,directionobs,len,epsilon2;
        eps2_boundary_constrain = eps2_boundary_constrain,
        eps2_div_constrain = eps2_div_constrain,
        eps2_Coriolis_constrain = eps2_Coriolis_constrain,
        f = f,
        g = g,
        residual = residual,
        lenη = lenη,
        ratio = ratio
    )


    residual4 = fill(NaN,(imax,jmax,length(ntimes),nsites))
    residual4[valid] = residual
    #res[:,:,ncenter,:] = residual4[:,:,(end+1) ÷ 2,:]

    #return residual4[:,:,(end+1) ÷ 2,:]
    return uri_temp[:,:,(end+1) ÷ 2],vri_temp[:,:,(end+1) ÷ 2],ηri_temp[:,:,(end+1) ÷ 2],residual4[:,:,(end+1) ÷ 2,:]
end

function cverr(
    xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
    lonr,latr,timerange,
    mask2d,h,
    len,lenη,eps2,
    eps2_boundary_constrain,
    eps2_div_constrain,
    eps2_Coriolis_constrain,
    g,ratio; u = [], v = [], η = [], selection = :cv)


    # time window
    Δn = 1

    # DEBUG only 1:3

    if selection == :all
        ncv = collect(1:size(flagcv_all,3))
    elseif selection == :cv
        ncv = findall(sum(flagcv_all,dims = [1,2,4])[:] .> 0)
    elseif selection == :debug
        @show "only 2"

        ncv = findall(sum(flagcv_all,dims = [1,2,4])[:] .> 0)[1:2]
    else
        error("unknown selection")
    end

    #ncv = findall(sum(flagcv_all,[1 2 4])[:] .> 0)

    fun(ncenter) = cv_rec(ncenter,Δn,
                          xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,sitenames,
                          lonr,latr,timerange,
                          mask2d,h,
                          len,lenη,eps2,
                          eps2_boundary_constrain,
                          eps2_div_constrain,
                          eps2_Coriolis_constrain,
                          g,
                          ratio
                          )

    uri = Array{Float64}(undef,(length(lonr),length(latr),length(timerange)))
    vri = Array{Float64}(undef,(length(lonr),length(latr),length(timerange)))
    ηri = Array{Float64}(undef,(length(lonr),length(latr),length(timerange)))
    res = Array{Float64}(undef,size(robs_all))

    uri[:] .= NaN
    vri[:] .= NaN
    ηri[:] .= NaN
    res[:] .= NaN

    #@show fun(50)
    #@show ncv
    output = pmap(fun,ncv)
    #output = map(fun,ncv)
    #@show size(output),size(ncv)

    uri[:,:,ncv] = cat([o[1] for o in output]...; dims = 4)
    vri[:,:,ncv] = cat([o[2] for o in output]...; dims = 4)
    ηri[:,:,ncv] = cat([o[3] for o in output]...; dims = 4)
    res[:,:,ncv,:] = permutedims(cat([o[4] for o in output]...; dims = 4),[1,2,4,3])

    #@show size(res),size(flagcv_all),size(output)
    # 0.0652580579558992 without Coriolis and no div
    # 0.05105108269989013 with f and g, no div

    # over all RMS
    forcv_and_sea = flagcv_all .& (.!isnan.(res))
    cv_err = sqrt(mean(res[forcv_and_sea].^2))
    status = (len...,lenη...,eps2,
              eps2_boundary_constrain,
              eps2_div_constrain,
              eps2_Coriolis_constrain,
              g,ratio,
              cv_err)

    @show status

    open("temp-output$(get(ENV,"SLURM_JOB_NAME",""))$(get(ENV,"SLURM_JOBID","")).txt","a") do f
        print(f,join(status,'\t'),"\n")
    end

    #info("over all cross-validation central-time error is $(cv_err)")

    if !isempty(u)
        u[:] = uri
    end

    if !isempty(v)
        v[:] = vri
    end

    if !isempty(η)
        η[:] = ηri
    end

    GC.gc()
    return cv_err
end
