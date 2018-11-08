struct Gaussian{DIM, T <:AbstractFloat}
    σ::T
    x0::Vec{DIM, T}
end
Gaussian(σ, crystal::Crystal) = Gaussian(σ, cellcenter(crystal))
(x::Gaussian{DIM})(r::Vec{DIM}) where DIM =  1/(2π * x.σ^2)^(DIM/2) * exp(-norm(r - x.x0)^2 / (2x.σ^2))
