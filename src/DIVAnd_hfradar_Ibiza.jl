
println("Ibiza_case")

Ibiza_case = (
    # string identifier of the sites
    sites = ["FORM","GALF"],
    # contains the file radials.nc
    basedir = joinpath(ENV["HOME"],"tmp/HFRadar-Ibiza/Radials/"),
    # bathymetry
    bathname = joinpath(ENV["HOME"],"tmp","HFRadar-Ibiza","bathymetry.mat")
    # directory for the figures
    figdir = joinpath(ENV["HOME"],"Doc/DIVAnd_hfradar/Fig"),
#    imax = 100,
#    jmax = 100,
    lonrange = [0.3,1.45],
    latrange = [38.3,39.4],
    # reduction factor in resolution of the bathymetry
    red = 16,
    # maximum depth of surface currents
    hmax = 50,
)


mycase = Ibiza_case
