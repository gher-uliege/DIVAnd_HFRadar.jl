

function hfradar_plot(
    xi,yi,datetime,uri,vri,
    xobs,yobs,robs,directionobs,sitenames,siteorigin;
    cmap = "jet",
    fig_suptitle = ""
)



    fig = figure(figsize=(10,6))
    #subplots_adjust(hspace=0., wspace=0.4)
#    subplots_adjust(hspace=0.1, wspace=0.4)
    subplots_adjust(wspace=0.35)
    fig.suptitle("$fig_suptitle - $(datetime)")

    #subplot(2,2,1)
    subplot(1,2,1)

    color_quiver_analysis = "#3A81CC"
    color_quiver_legend = "#C44E52"


    scale = 7
    scale = 5
    quiver(xi,yi,uri,vri,color=color_quiver_analysis,scale = scale)

    title("(a) Total currents")
    plotcoastline()
    quiver([1.3],[38.95],[0.5],[0.],color=color_quiver_legend,scale = scale)

    xlim(lonr[[1,end]])
    ylim(latr[[1,end]])
    set_aspect_ratio()
    gca().set_xticks([0.5,1])

    #quiver(xobs,yobs,uobs,vobs,color="r",scale = scale);

    cl = extrema(robs[isfinite.(robs)])
    cldiff = [-0.15,0.15]
@show cl
ax1 = nothing
ax_diff = nothing
for n = 1:length(sitenames)
#    global ax1
#    global ax_diff
        site = sitenames[n]
        sel = (:,:,1,n)
        @show site

        longitude,latitude = siteorigin[site]

        bearingi = GeoMapping.azimuth.(latitude,longitude,yi,xi)
        directioni = bearingi .+ 180

        ri = uri .* sin.(directioni * pi/180) + vri .* cos.(directioni * pi/180)
        itp = extrapolate(interpolate((xi[:,1],yi[1,:]),ri, Gridded(Linear())),NaN)
        ri_at_obs = itp.(xobs[:,:,n],yobs[:,:,n])

        #subplot(2,2,2+n)
        ax_velocity = subplot(2,4,2+n)
        pcolor(xi,yi,ri; cmap=cmap)

        plotcoastline()
        #scatter(xobsn,yobsn,20,robsn,edgecolor="none")
        ax1 = scatter(xobs[:,:,n],yobs[:,:,n],4,robs[:,:,n]; cmap=cmap)
        clim(cl...)
        set_aspect_ratio()
        #colorbar()
        xlim(lonr[[1,end]])
        ylim(latr[[1,end]])
        title("($('a' + n)) Rad. curr. - $(site)")

        subplot(2,4,2+n+4)
        plotcoastline()
        #scatter(xobsn,yobsn,20,robsn,edgecolor="none")
        ax_diff = scatter(xobs[:,:,n],yobs[:,:,n],10,ri_at_obs - robs[:,:,n]; cmap="bwr")
        clim(cldiff...)
        set_aspect_ratio()
        #colorbar()
        xlim(lonr[[1,end]])
        ylim(latr[[1,end]])
        title("($('a' + n + 2)) Diff. - $(site)")
    end

    cb_w,cb_h = (0.015, ax1.axes.get_position().height)
    cb_r = 0.95
    cb_r = 0.92
    cbaxes = fig.add_axes([cb_r, ax1.axes.get_position().y0, cb_w, cb_h])
    colorbar(ax1,cax = cbaxes,orientation="vertical")

    cbaxes = fig.add_axes([cb_r, ax_diff.axes.get_position().y0, cb_w, cb_h])
    colorbar(ax_diff,cax = cbaxes,orientation="vertical")
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

