using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
using Documenter, DIVAnd_HFRadar

makedocs(modules = [DIVAnd_HFRadar], sitename = "DIVAnd_HFRadar.jl")

if CI
    deploydocs(repo = "github.com/gher-ulg/DIVAnd_HFRadar.jl.git",
               target = "build")
end
