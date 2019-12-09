""""
Stagger from a `v` to a `rho` location in an Arakawa C grid
"""
function stagger_v2r(v)

    sz = [size(v)...]
    sz[2] = sz[2]+1
    r = zeros((sz...,))

    n = size(v,2)+1

    r[:,2:n-1,:,:] = (v[:,1:n-2,:,:]+v[:,2:n-1,:,:])/2
    r[:,1,:,:] = v[:,1,:,:]
    r[:,n,:,:] = v[:,n-1,:,:]

    return r
end
