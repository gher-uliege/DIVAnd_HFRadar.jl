#!/usr/bin/env julia
using Glob
using Base.Markdown

include("DIVAnd_hfradar_load_optim.jl")


ref = 0.; 
postfix = ARGS[1]

cases = ["2D","2D_bc","2D_div","3D","3D_Coriolis","3D_Coriolis_geo","3D_Coriolis_geo_gp"]
rows = Any[Any["Case","RMS","Skill","Optimal parameter(s)"]]; 

for c in cases

    RMS,param = load_optim_param(c,postfix)

    if length(param) > 0
        
        op = ""
        pp = ["eps2"]
        if c == "2D"
            ref = RMS
        elseif c == "2D_boundary_constrain"
            pp = ["eps2","eps2_boundary_constrain"]
        elseif c == "2D_div_constrain"
            pp = ["eps2","eps2_div_constrain"]
        elseif c == "3D"
            pp = ["eps2","lent"]
        elseif c == "3D_Coriolis_constrain"
            pp = ["eps2","eps2_Coriolis_constrain"]
        elseif c == "3D_Coriolis_geo" || c == "3D_Coriolis_geo_gp"
            pp = ["eps2","eps2_Coriolis_constrain","ratio"]
        end 
        
        oplist = join(["$(p)=$(@sprintf("%.4g",param[p]))" for p in pp],", ")
        push!(rows,[c,@sprintf("%.4f",RMS), @sprintf("%.3f",1 - RMS^2/ref^2), oplist])

        @show c,param
    end
end
    
print(Base.Markdown.plain(Base.Markdown.Table(rows,[:l,:r,:r,:l])))