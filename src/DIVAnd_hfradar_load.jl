using Compat: findall

const sites = ["FORM","GALF"]
const basedir   = joinpath(ENV["HOME"],"tmp/HFRadar-Ibiza/Radials/")


const imax = 100;
const jmax = 100;
function daNA(tmp)
    return nomissing(tmp,NaN)
end

#peshape(xobs_all,(size(xobs_all,1),size(xobs_all,2),1,size(xobs_all,3)))

const filename = joinpath(basedir,"radials.nc")

ds = Dataset(filename)
const timerange = nomissing(ds["time"][:])
const robs_all = daNA(ds["radial"][:])
const xobs_all = daNA(ds["lon"][:])
const yobs_all = daNA(ds["lat"][:])
const directionobs_all = daNA(ds["direction"][:])
const flagcv_all = Bool.(nomissing(ds["flagcv"][:]))

const sitenames = nomissing(ds["site"][:])
const sites_lon = nomissing(ds["site_lon"][:])
const sites_lat = nomissing(ds["site_lat"][:])


const figdir = joinpath(ENV["HOME"],"Doc/DIVAnd_hfradar/Fig")

if false
    clf();
    plot(time_all,sum(.!isnan.(robs_all),[1 2 4])[:]; label="number of radials")
    plot(time_all,sum(flagcv_all,[1 2 4])[:]; label="number of cross-validation points");
    title("Temporal statistics")
    legend()
    savefig(joinpath(figdir,"nb_cv_points.eps"))
    savefig(joinpath(figdir,"nb_cv_points.svg"))
end

const siteorigin = Dict(sitenames[i] => (sites_lon[i],sites_lat[i]) for i = 1:length(sitenames))

#filenames = [joinpath(basedir,"$(site)/$(Dates.format(t,"yyyy_mm"))/RDLm_$(site)_$(Dates.format(t,"yyyy_mm_dd_HHMM")).ruv") for t in timerange, site in sites]
#xobs_all,yobs_all,robs_all,directionobs_all,flagcv_all,metadata = loadall(filenames,imax,jmax)
# convert from cm/s to m/s
#robs_all = robs_all/100
# robs_all = robs_all/100

# sitenames = Vector{String}(length(sites))
# siteorigin = Dict{String,Tuple{Float64,Float64}}()

# for isite = 1:length(sites)
#     sitenames[isite] = strip(replace(metadata[isite]["Site"],"\"\"",""))
#     latitude,longitude = [parse(Float64,c) for c in split(metadata[isite]["Origin"])]
#     siteorigin[sitenames[isite]] = (longitude,latitude)
# end








# lon/lat range

const lonrange = [0.3,1.45]
const latrange = [38.3,39.4]


#mxi,myi,mask = load_mask(bath_name,true,lonr,latr,0)
bathname = joinpath(basedir,"../bathymetry.mat")
origh = MAT.matread(bathname)["bat"]
origlon = MAT.matread(bathname)["lon"]
origlat = MAT.matread(bathname)["lat"]

j = findall((latrange[1] .< origlat) .& (origlat .< latrange[end]))
i = findall((lonrange[1] .< origlon) .& (origlon .< lonrange[end]))

red = 16

const lonr = origlon[i[1]:red:i[end]]
const latr = origlat[j[1]:red:j[end]]

htmp = origh[i[1]:red:i[end],j[1]:red:j[end]]
const mask2d = htmp .< 0
htmp[.!mask2d] .= 0
htmp = -htmp
#const h = min.(htmp,75.)
const h = min.(htmp,50.)
const htot = htmp



function plotcoastline()
    j = findall((latr[1] .< origlat) .& (origlat .< latr[end]))
    i = findall((lonr[1] .< origlon) .& (origlon .< lonr[end]))
    contourf(origlon[i],origlat[j],origh[i,j]' .< 0,levels=[0,.5], cmap = "gray")
end


# rotation of earth (siderial day)
const Ω = 7.2921e-5 # rad/s
const f = 2*Ω*sin(mean(latr) * π /180)

#ρ_1 = PhysOcean.density0(38.2,13) 
#ρ_0 = PhysOcean.density0(37.7,20)
ρ_0 = 1026.6;
ρ_1 = 1028.7;

const g_barotropic = 9.81
const g_baroclinic = g_barotropic * 2 * (ρ_1 - ρ_0) / (ρ_1 + ρ_0)


#g = 0.


