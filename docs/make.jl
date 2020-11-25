using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
using Documenter, DIVAnd_hfradar

makedocs(modules = [DIVAnd_hfradar], sitename = "DIVAnd_hfradar.jl")

if CI
    deploydocs(repo = "github.com/gher-ulg/DIVAnd_hfradar.jl.git",
               target = "build")
end
