
Eself(ρ::Gaussian) = 1/(2*√π * ρ.σ)
function Eself(crystal::Crystal, ρ::Gaussian{T}) where T
    E = zero(T)
    for at in eachatom(crystal)
        E += at[:charge] * Eself(ρ)
    end
    return E
end

EHartree(ϕ, ρ, grid::PWGrid) = sum(ϕ .* ρ) * dΩ(grid)/2

function Eewald(ρ, grid, Eself)
    ϕ = solve_poisson(ρ, grid) #potential from ρ
    return EHartree(ϕ, ρ, grid) - Eself
end
