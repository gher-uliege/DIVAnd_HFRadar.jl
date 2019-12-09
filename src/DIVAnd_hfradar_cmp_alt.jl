using NCDatasets
using OceanPlot
using PyPlot

fnames = OceanPlot.listfiles(expanduser("~/Data/Med/Alt/"); extension = "nc")
#fnames = OceanPlot.listfiles(expanduser("~/Data/Med/Alt2/al/"); extension = "nc")
#fnames = OceanPlot.listfiles(expanduser("~/Data/Med/Alt2/all/"); extension = "nc")


radfname = "/media/abarth/9982a2e8-599f-4b37-884f-a59aa1c0de80/Alex/nic4/tmp/HFRadar-Ibiza/3D_Coriolis_geo.nc"
radfname = expanduser("~/tmp/HFRadar-Ibiza/3D_Coriolis_geo.nc")
radfname = expanduser("~/tmp/HFRadar-Ibiza/rg/3D_Coriolis_geo.nc")
radfname = expanduser("~/tmp/HFRadar-Ibiza/rg/test.nc")

ds = Dataset(radfname);

i = 10:35;
j = :;

x = nomissing(ds["lon"][i]);
y = nomissing(ds["lat"][j]);
deta = ds["eta"][i,j,:];



cl = (-0.02,0.06)
cmap = "jet"

#for sat in ["al","c2","h2","j2"]

#    figure()
#fnames = OceanPlot.listfiles(expanduser("~/Data/Med/Alt2/$sat/"); extension = "nc")

mean_eta = mean(deta,3)[:,:,1]


pcolor(x,y,mean(deta,3)[:,:,1]'; vmin = cl[1], vmax = cl[2], cmap= cmap)
for fname in fnames
    @show fname
    ds = Dataset(fname)
    lon = nomissing(ds["longitude"][:])
    lat = nomissing(ds["latitude"][:]);

    lon = mod.(lon + 180,360) - 180;

    sl = ds["adt_filtered"][:];
    #sl = ds["adt_unfiltered"][:];
    scatter(lon,lat,20,sl; vmin = cl[1], vmax = cl[2], cmap= cmap)
    close(ds)
end


colorbar()
xlim(rg(x)...)
ylim(rg(y)...)
#    title(sat)
#end
