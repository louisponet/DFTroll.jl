module DFTroll
    using LinearAlgebra, Base.Cartesian
    using StaticaArrays
    using FFTW, Unitful, Crystals
    using Parameters
    # Unitful.register(@__MODULE__)
    include("constants.jl")
    # include("types.jl")
        # const UH = UnitfulHartree
    # export UnitfulHartree, Unitful, @u_str
    # export a₀

    include("utils.jl")
    include("grids.jl")
    include("distributions.jl")
    include("energy.jl")
    include("calculation.jl")
    export Calculation, Gaussian, Crystal
    export pwcalculation, density_integral

    const localunits = Unitful.basefactors
    function __init__()
        merge!(Unitful.basefactors, localunits)
        Unitful.register(DFTroll)
    end
    export a₀
end # module
