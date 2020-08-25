

function hfradar_plot(
    xi,yi,datetime,uri,vri,
    xobs,yobs,robs,directionobs,sitenames,siteorigin;
    cmap = "jet"
)

    figure(figsize=(9,4))
    subplots_adjust(hspace=0., wspace=0.4)

    #subplot(2,2,1)
    subplot(1,3,1)

    color_quiver_analysis = "#3A81CC"
    color_quiver_legend = "#C44E52"


    scale = 7
    quiver(xi,yi,uri,vri,color=color_quiver_analysis,scale = scale)

    title("$(datetime)")
    plotcoastline()
    quiver([1.3],[38.95],[0.5],[0.],color=color_quiver_legend,scale = scale)

    xlim(lonr[[1,end]])
    ylim(latr[[1,end]])
    set_aspect_ratio()

    #quiver(xobs,yobs,uobs,vobs,color="r",scale = scale);

    cl = extrema(robs[isfinite.(robs)])
    @show cl
    for n = 1:length(sitenames)
        site = sitenames[n]
        sel = (:,:,1,n)
        @show site

        longitude,latitude = siteorigin[site]

        bearingi = GeoMapping.azimuth.(latitude,longitude,yi,xi)
        directioni = bearingi .+ 180

        ri = uri .* sin.(directioni * pi/180) + vri .* cos.(directioni * pi/180)

        #subplot(2,2,2+n)
        subplot(1,3,1+n)
        pcolor(xi,yi,ri; cmap=cmap)

        plotcoastline()
        #scatter(xobsn,yobsn,20,robsn,edgecolor="none")
        ax1 = scatter(xobs[:,:,n],yobs[:,:,n],10,robs[:,:,n]; cmap=cmap)
        clim(cl...)
        set_aspect_ratio()
        #colorbar()
        xlim(lonr[[1,end]])
        ylim(latr[[1,end]])
        title("Radial currents $(site)")
    end

    fig = gcf()
    cbaxes = fig.add_axes([0.4, 0.1, 0.5, 0.03])
    colorbar(cax = cbaxes,orientation="horizontal")

#    savefig("DIVAnd_hfradar_2D.png"; dpi=200)
#    savefig("DIVAnd_hfradar_2D.svg"; dpi=200)

end



function plotradar(x,y,robs,flagcv,sitenames,time)
    clf()
    for j = 1:size(robs,3)
        @show size(robs),size(x)
        subplot(2,2,j); scatter(x[:,:,j],y[:,:,j],10,robs[:,:,j]; cmap="jet"); title("$(sitenames[j]) $(time)")
        plotcoastline()
        set_aspect_ratio()

        xlim(extrema(xobs_all[isfinite.(xobs_all)]))
        ylim(extrema(yobs_all[isfinite.(yobs_all)]))

        subplot(2,2,2+j); scatter(x[:,:,j],y[:,:,j],10,flagcv[:,:,j]; cmap="jet"); title("flag $(sitenames[j])")
        plotcoastline()
        xlim(extrema(xobs_all[isfinite.(xobs_all)]))
        ylim(extrema(yobs_all[isfinite.(yobs_all)]))


        set_aspect_ratio()
    end

    savefig("Fig/plotradar_$(time).png")

end
