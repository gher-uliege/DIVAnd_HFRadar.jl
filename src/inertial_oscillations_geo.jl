""" 
u size imax-1,jmax and v of size imax,jmax-1; velocities on a Arakawa C grid
(time invariant)
"""
struct CGrid{TA <: AbstractArray{<:Number,2},BA <: AbstractArray{Bool}}
    mask::BA
    mask_u::BA
    mask_v::BA

    h::TA
    h_u::TA
    h_v::TA
    
    pm::TA
    pm_u::TA
    pm_v::TA

    pn::TA
    pn_u::TA
    pn_v::TA
    
end

# model configuration
struct Config{TA,BA,T}
    grid::CGrid{TA,BA}
    f::T
    g::T
end

function CGrid(mask,h,pmn)
    # staggering

    mask_u, mask_v, mask_psi = stagger_mask(mask)

    h_u = stagger_r2u(h)
    h_v = stagger_r2v(h)

    pmn_u = [stagger_r2u(pm) for pm in pmn]    
    pmn_v = [stagger_r2v(pm) for pm in pmn]

    return CGrid(mask,mask_u,mask_v,
                 h,h_u,h_v,
                 pmn[1],pmn_u[1],pmn_v[1],
                 pmn[2],pmn_v[2],pmn_v[2])    
end

# function Config_(grid,Δt,f,g)
#     return Config(grid,Δt,f,g)
# end

function intertial_oscillation_geo(config,Δt,η,u,v,δη,δu,δv)

    δηn = zeros(size(δv))
    δun = zeros(size(δu))
    δvn = zeros(size(δv))

    δηn = δη
    Δtf = config.f * Δt
    Δtg = config.g * Δt
    
    # update δu
    for i = 1:size(δu,1)
        for j = 1:size(δu,2)   
            jc = min(j,size(δv,2))
            jm = max(j-1,1)
            δvp = (δv[i,jc] + δv[i+1,jc] + δv[i,jm] + δv[i+1,jm])/4

            δun[i,j] = cos(Δtf) * δu[i,j] + sin(Δtf) * δvp -
                Δtg * config.grid.pm_u[i,j] * (δη[i+1,j] - δη[i,j])
        end
    end


    # update δv
    for i = 1:size(δv,1)
        for j = 1:size(δv,2)
            ic = min(i,size(δu,1))
            im = max(i-1,1)
            δup = (δu[ic,j] + δu[im,j] + δu[ic,j+1] + δu[im,j+1])/4

            δvn[i,j] = cos(Δtf) * δv[i,j] - sin(Δtf) * δup -
                Δtg * config.grid.pn_v[i,j] * (δη[i,j+1] - δη[i,j])

        end
    end

    return  (δηn,δun,δvn)
end

function intertial_oscillation_geo_adj(config,Δt,η,u,v,δηn,δun,δvn)
    δη = zeros(size(δηn))
    δu = zeros(size(δun))
    δv = zeros(size(δvn))

    δη = δηn
    Δtf = config.f * Δt
    Δtg = config.g * Δt
    
    # update δv
    for i = 1:size(v,1)
        for j = 1:size(v,2)
            #vn = v + Δtf *up
            δv[i,j] = δv[i,j] + cos(Δtf) * δvn[i,j]
            δup = -sin(Δtf) * δvn[i,j]

            ic = min(i,size(u,1))
            im = max(i-1,1)

            δu[ic,j] = δu[ic,j] + δup/4
            δu[im,j] = δu[im,j] + δup/4
            δu[ic,j+1] = δu[ic,j+1] + δup/4
            δu[im,j+1] = δu[im,j+1] + δup/4

            tmp = Δtg * config.grid.pn_v[i,j]
            δη[i,j+1] -=  tmp * δvn[i,j]
            δη[i,j] += tmp * δvn[i,j]

        end
    end


    # update δu
    for i = 1:size(u,1)
        for j = 1:size(u,2)            
            #un = u - Δtf * vp
            δu[i,j] = δu[i,j] + cos(Δtf) * δun[i,j]
            δvp = +sin(Δtf) * δun[i,j]
 
            jc = min(j,size(v,2))
            jm = max(j-1,1)

            δv[i,jc] = δv[i,jc] + δvp/4
            δv[i+1,jc] = δv[i+1,jc] + δvp/4
            δv[i,jm] = δv[i,jm] + δvp/4
            δv[i+1,jm] = δv[i+1,jm] + δvp/4

            tmp = Δtg * config.grid.pm_u[i,j]
            δη[i+1,j] -=  tmp * δun[i,j]
            δη[i,j] +=  tmp * δun[i,j]

        end
    end
    return  (δη,δu,δv)
end


function misfit_intertial_oscillation_geo(sv,config,t,x,δx)    
    η,u,v = DIVAnd.unpack(sv,x)
    δη,δu,δv = DIVAnd.unpack(sv,δx)
    # leave first instance equal to zero
    δdiffη = zeros(size(η))
    δdiffu = zeros(size(u))
    δdiffv = zeros(size(v))

    for n = 2:size(u,3)
        #Δtf = config.f * (t[n]-t[n-1])
        Δt = t[n]-t[n-1]
        δηp,δup,δvp = intertial_oscillation_geo(config,Δt,
                                            η[:,:,n-1],u[:,:,n-1],v[:,:,n-1],
                                            δη[:,:,n-1],δu[:,:,n-1],δv[:,:,n-1])
        δdiffη[:,:,n] = δηp - δη[:,:,n]
        δdiffu[:,:,n] = δup - δu[:,:,n]
        δdiffv[:,:,n] = δvp - δv[:,:,n]
    end

    return DIVAnd.pack(sv,(δdiffη,δdiffu,δdiffv))[:,1]
end


function misfit_intertial_oscillation_geo_adj(sv,config,t,x,δdiffx)
    η,u,v = DIVAnd.unpack(sv,x)
    δdiffη,δdiffu,δdiffv = DIVAnd.unpack(sv,δdiffx)

    δη = zeros(size(δdiffη))
    δu = zeros(size(δdiffu))
    δv = zeros(size(δdiffv))

    for n = size(δdiffu,3):-1:2
        Δt = t[n]-t[n-1]
        δη[:,:,n] = δη[:,:,n] - δdiffη[:,:,n]
        δu[:,:,n] = δu[:,:,n] - δdiffu[:,:,n]
        δv[:,:,n] = δv[:,:,n] - δdiffv[:,:,n]
        
        δη[:,:,n-1],δu[:,:,n-1],δv[:,:,n-1] = intertial_oscillation_geo_adj(config,
            Δt,
            η[:,:,n-1],u[:,:,n-1],v[:,:,n-1],
            δdiffη[:,:,n],δdiffu[:,:,n],δdiffv[:,:,n])
    end

    return DIVAnd.pack(sv,(δη,δu,δv))[:,1]
end

