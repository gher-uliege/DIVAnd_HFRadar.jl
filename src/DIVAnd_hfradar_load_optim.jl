function load_optim_param(c,postfix)
    files = sort(glob("$(c)-$(postfix)/*txt"))

    @show files
    if length(files) > 0
        data = readdlm(files[end])
        #@show data[findmin(data[:,end])[2],:]
        n = findmin(data[:,end])[2]
        RMS = data[n,end]

        d = data[n,:]

        return RMS,Dict("lenxy" => d[1],"lent" => d[3], "lenetat" => d[6], "eps2" => d[7], "eps2_boundary_constrain" => d[8],
                    "eps2_div_constrain" => d[9], "eps2_Coriolis_constrain" => d[10], "g" => d[11], "ratio" => d[12])
    end

    return NaN,Dict()
end
