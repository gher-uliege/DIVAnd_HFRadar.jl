""" 
u size imax-1,jmax and v of size imax,jmax-1; velocities on a Arakawa C grid
"""
function intertial_oscillation(Δtf,u,v,δu,δv)

    δun = zeros(size(δu))
    δvn = zeros(size(δv))

    # update δu
    for i = 1:size(δu,1)
        for j = 1:size(δu,2)   
            jc = min(j,size(δv,2))
            jm = max(j-1,1)
            δvp = (δv[i,jc] + δv[i+1,jc] + δv[i,jm] + δv[i+1,jm])/4

            δun[i,j] = cos(Δtf) * δu[i,j] + sin(Δtf) * δvp
        end
    end


    # update δv
    for i = 1:size(δv,1)
        for j = 1:size(δv,2)
            ic = min(i,size(δu,1))
            im = max(i-1,1)
            δup = (δu[ic,j] + δu[im,j] + δu[ic,j+1] + δu[im,j+1])/4

            δvn[i,j] = cos(Δtf) * δv[i,j] - sin(Δtf) * δup
        end
    end

    return  (δun,δvn)
end

function intertial_oscillation_adj(Δtf,u,v,δun,δvn)
    δu = zeros(size(δun))
    δv = zeros(size(δvn))
    
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
        end
    end
    return  (δu,δv)
end


function misfit_intertial_oscillation(sv,f,t,x,δx)    
    u,v = DIVAnd.unpack(sv,x)
    δu,δv = DIVAnd.unpack(sv,δx)
    # leave first instance equal to zero
    δdiffu = zeros(size(u))
    δdiffv = zeros(size(v))

    for n = 2:size(u,3)
        Δtf = f * (t[n]-t[n-1])
        δup,δvp = intertial_oscillation(Δtf,u[:,:,n-1],v[:,:,n-1],δu[:,:,n-1],δv[:,:,n-1])
        δdiffu[:,:,n] = δup - δu[:,:,n]
        δdiffv[:,:,n] = δvp - δv[:,:,n]
    end

    return DIVAnd.pack(sv,(δdiffu,δdiffv))[:,1]
end


function misfit_intertial_oscillation_adj(sv,f,t,x,δdiffx)
    u,v = DIVAnd.unpack(sv,x)
    δdiffu,δdiffv = DIVAnd.unpack(sv,δdiffx)

    δu = zeros(size(δdiffu))
    δv = zeros(size(δdiffv))

    for n = size(δdiffu,3):-1:2
        Δtf = f * (t[n]-t[n-1])
        δu[:,:,n] = δu[:,:,n] - δdiffu[:,:,n]
        δv[:,:,n] = δv[:,:,n] - δdiffv[:,:,n]
        
        δu[:,:,n-1],δv[:,:,n-1] = intertial_oscillation_adj(Δtf,u[:,:,n-1],v[:,:,n-1],δdiffu[:,:,n],δdiffv[:,:,n])
    end

    return DIVAnd.pack(sv,(δu,δv))[:,1]
end

