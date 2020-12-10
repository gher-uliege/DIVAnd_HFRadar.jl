using NCDatasets

function DIVAnd_HFRadar_save(fname,lonr,latr,timerange,uri,vri,ηri)

    function MArray(data)
        tmp = Array{Union{eltype(data),Missing}}(data)
        tmp[isnan.(tmp)] .= missing
        return tmp
    end

    fv = -99999.

    Dataset(fname,"c") do nc
        defDim(nc,"lon",size(uri,1))
        defDim(nc,"lat",size(uri,2))
        defDim(nc,"time",size(uri,3))

        nclon = defVar(nc,"lon",Float64,("lon",))
        nclat = defVar(nc,"lat",Float64,("lat",))
        nctime = defVar(nc,"time",Float64,("time",), attrib = Dict(
            "units" => "days since 1900-01-01 00:00:00"
        ))

        ncu = defVar(nc,"u",Float64,("lon","lat","time"); fillvalue = fv)
        ncv = defVar(nc,"v",Float64,("lon","lat","time"); fillvalue = fv)
        nceta = defVar(nc,"eta",Float64,("lon","lat","time"); fillvalue = fv)

        nclon[:] = lonr
        nclat[:] = latr
        nctime[:] = timerange

        ncu[:] = MArray(uri)
        ncv[:] = MArray(vri)
        nceta[:] = MArray(ηri)
    end

    return nothing
end

export DIVAnd_HFRadar_save
