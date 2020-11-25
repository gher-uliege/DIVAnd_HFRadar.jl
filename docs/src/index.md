# DIVAnd HF Radar


The package `DIVAnd_hfradar` allow to interpolate surface current data on a regular grid.
The primary use-case is for radial current measurements for high-frequency radars. But it can also be applied to any other
current data (like ADCPs or drifters).


## Tutorial

To run these examples you need install Julia and the packages `DIVAnd`, `DIVAnd_hfradar` and `PyPlot` which can be installed by
these julia commands:

```julia
using Pkg
Pkg.add("DIVAnd")
Pkg.add(url="https://github.com/gher-ulg/DIVAnd_hfradar.jl", rev="master")
Pkg.add("PyPlot")
```

In Linux, you must also install the python package `matplotlib`. Under Debian/Ubuntu, you can do this via the shell command:

```bash
sudo apt install python3-matplotlib
```

### Data constraint

In this example, we are setting up an idealized domain from spanning -1 to 1 with 10x11 grid points.
The gray area on the right is a coastal wall.

```@example 1
using DIVAnd_hfradar: DIVAndrun_hfradar
using DIVAnd
using PyPlot

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

# helper function to plot results
function plotres(uri,vri)
    clf()
    quiver(xi,yi,uri,vri, scale = 10)
    α = directionobs*pi/180
    quiver(xobs,yobs,robs .* sin.(α), robs .* cos.(α),color = "r",scale = 10)
    contourf(xi,yi,mask,levels = [0,.5],cmap = "gray")
end

uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn),(xi,yi),(xobs,yobs),robs,directionobs,len,epsilon2)
plotres(uri,vri)
title("Data constraint")
savefig("currents1.png"); clf(); nothing # hide
```

![](currents1.png)


### Boundary condition

```@example 1
uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn),(xi,yi),(xobs,yobs),robs,directionobs,len,epsilon2,
    eps2_boundary_constrain = 0.0001,
)
plotres(uri,vri)
title("Data constraint and boundary condition")
savefig("currents2.png"); clf(); nothing # hide
```

![](currents2.png)

### Boundary condition and divergence constraint

```@example 1
uri,vri = DIVAndrun_hfradar(
    mask,h,(pm,pn),(xi,yi),(xobs,yobs),robs,directionobs,len,epsilon2,
    eps2_boundary_constrain = 0.001,
    eps2_div_constrain = 0.001,
)
plotres(uri,vri)
title("Data constraint, boundary condition and divergence constraint")
savefig("currents3.png"); clf(); nothing # hide
```

![](currents3.png)

### 3D anaysis with time and the Coriolis force

```@example 1
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

α = directionobs*pi/180

quiver(xi[:,:,1],yi[:,:,1],uri[:,:,1],vri[:,:,1],
       color="#b2c4db",label="t-1 hour");
quiver(xi[:,:,2],yi[:,:,2],uri[:,:,2],vri[:,:,2],
       color="k",label="current time t");
quiver(xi[:,:,3],yi[:,:,3],uri[:,:,3],vri[:,:,3],
       color="#b2dbbd",label="t+1 hour");
quiver(xobs,yobs,robs .* sin.(α), robs .* cos.(α),
       color = "#e86966",scale = 10,label="rad. obs. at time t")
legend(loc="upper right")

savefig("currents_coriolis.png"); nothing # hide
```

![](currents_coriolis.png)


## Reference

```@autodocs
Modules = [DIVAnd_hfradar]
Order   = [:function, :type]
```

