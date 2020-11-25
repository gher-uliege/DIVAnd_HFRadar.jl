using DIVAnd_hfradar: DIVAndrun_hfradar
using DIVAnd

# size of the grid
sz = (10,11)

# depth (meters)
h = 50 * ones(sz)

# land-sea mask
# true is sea; false is land
mask = trues(sz)
mask[end,:] .= false

# 2D grid
xi,yi = DIVAnd.ndgrid(LinRange(-1,1,sz[1]),LinRange(-1,1,sz[2]))

# scale factor; inverse of the resolution
pm = ones(sz) / (xi[2,1]-xi[1,1]);
pn = ones(sz) / (yi[1,2]-yi[1,1]);

# radial observations
robs = [1.]

# direction of the observation (from North counted clockwise)
directionobs = [90.]

# position of the observation
xobs = [0.]
yobs = [0.]

# correlation length
len = (0.6,0.6)

# data constraint
epsilon2 = 0.001

uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn),(xi,yi),(xobs,yobs),robs,directionobs,len,epsilon2)


uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn),(xi,yi),(xobs,yobs),robs,directionobs,len,epsilon2,
    eps2_boundary_constrain = 0.0001,
)


uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn),(xi,yi),(xobs,yobs),robs,directionobs,len,epsilon2,
    eps2_boundary_constrain = 0.001,
    eps2_div_constrain = 0.001,
)



# 3D grid (longitude, latitude and time)
sz = (10,11,3)

# depth
h = 50 * ones(sz)

# only see points
mask = trues(sz)

# 3D grid (time is in seconds)
xi,yi,ti = DIVAnd.ndgrid(
    range(-1,stop=1,length=sz[1]),
    range(-1,stop=1,length=sz[2]),
    range(-3600,stop=3600,length=sz[3]))

# scale factor; inverse of the resolution
pm = ones(size(xi)) / (xi[2,1,1]-xi[1,1,1]);
pn = ones(size(xi)) / (yi[1,2,1]-yi[1,1,1]);
po = ones(size(xi)) / (ti[1,1,2]-ti[1,1,1]);

# tobs is the time of the observations
# other variable are as before
robs = [1.]
xobs = [0.]
yobs = [0.]
tobs = [0.]
len = (0.6,0.6,0.0)
epsilon2 = 0.1

uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn,po),(xi,yi,ti),(xobs,yobs,tobs),robs,directionobs,len,epsilon2;
    eps2_boundary_constrain = -1,
    eps2_div_constrain = -1,
    eps2_Coriolis_constrain = 1e-1,
    f = 1e-4,
)

@test maximum(uri) > 0.7
