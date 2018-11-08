
LinearAlgebra.normalize(arr::AbstractArray) = arr./= norm(arr)

mutable struct PWWavefunction{T, DIM}
    k      ::Vec{DIM, T}
    coeffs ::Array{Complex{T},  DIM}
end
PWWavefunction(k::Vec{DIM,T} where DIM, g::NTuple) where T =
    PWWavefunction(k, normalize(rand(Complex{T}, g)))

ρreal(wfc::PWWavefunction) = real.(ifft(wfc.coeffs)).^2
ρrecip(wfc::PWWavefunction) = real.(fft(wfc.ρreal))

function orthogonalize!(wfcs::Vector{<:PWWavefunction})
    gridsize = size(wfcs[1].coeffs)
    flattenedcoeffs = hcat([reshape(wfc.coeffs, :) for wfc in wfcs]...)
    orthogonalize!(flattenedcoeffs)
    for i = 1:size(flattenedcoeffs)[2]
        wfcs[i].coeffs .= reshape(flattenedcoeffs[:, i], gridsize)
    end
end

LinearAlgebra.dot(wfc1::PWWavefunction, wfc2::PWWavefunction) = dot(wfc1.coeffs, wfc2.coeffs)




#Sets up all ion related and initializes the wavefunctions
struct Calculation{T, DIM}
    crystal      ::Crystal{T}
    grid         ::PWGrid{T, DIM}
    wavefunctions::Array{Vector{PWWavefunction{T, DIM}}, DIM}
    ρr           ::Array{T, DIM}
    Eewald       ::T
end
function pwcalculation(crystal::Crystal{T}, Ggridsize::NTuple{3}, kgridsize::NTuple{3}, ionic_ρdistribution, nbands, ) where  T
    pwgrid        = PWGrid(crystal.cell, Ggridsize)
    allwfcs =[[PWWavefunction(Vec{3, T}(k1, k2, k3), Ggridsize) for b=1:nbands] for k1=1:kgridsize[1],k2=1:kgridsize[2],k3=1:kgridsize[3]]

    for wfcs in allwfcs
        orthogonalize!(wfcs)
    end

    println(sum(sum(sum([ρreal.(wfcs) for wfcs in allwfcs]))))

    Sf = structure_factor(crystal, pwgrid)
    ρt = ionic_ρdistribution.(pwgrid.r) .+ 0.0im

    ρ = real.(ifft!(fft!(ρt) .* Sf))
    Ee = Eewald(ρ, pwgrid, Eself(crystal, ionic_ρdistribution))
    Calculation(crystal, pwgrid, allwfcs, ρ, Ee)
end

Ω(calc::Calculation) = Ω(calc.grid)
dΩ(calc::Calculation) = dΩ(calc.grid)
density_integral(calc::Calculation) = sum(calc.ρr) *  dΩ(calc)
