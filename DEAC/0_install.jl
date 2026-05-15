try
    using SmoQyDEAC
    using CairoMakie
catch
    using Pkg
    Pkg.add("SmoQyDEAC")
    Pkg.add("CairoMakie")
end