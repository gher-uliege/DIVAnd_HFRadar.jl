"""
Stagger from a `u` to a `rho` location in an Arakawa C grid
"""
function stagger_u2r(u)

    sz = [size(u)...]
    sz[1] = sz[1]+1
    r = zeros((sz...,))

    m = size(u,1)+1

    r[2:m-1,:,:,:] = (u[1:m-2,:,:,:]+u[2:m-1,:,:,:])/2
    r[1,:,:,:] = u[1,:,:,:]
    r[m,:,:,:] = u[m-1,:,:,:]

    return r
end
