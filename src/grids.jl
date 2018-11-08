
abstract type AbstractGrid{T, DIM} end

struct PWGrid{T, DIM} <: AbstractGrid{T, DIM}
    r ::VecArray{T, DIM} #Real space grid
    G ::VecArray{T, DIM} #Reciprocal space grid TODO: might be nice to have it centered around mid index?
    G2::Array{T, DIM} #absolute squared reciprocal space grid
    Ω ::T #Volume of the primitive cell
end
@generated function PWGrid(cell::Matrix{<:Quantity{T}}, gridsize::NTuple{N, Integer},) where {T, N}
    dimrange  = :(d -> 1:gridsize[d])
    preloop   = :(d -> (vec_d = Vec{N}(cell[:, d] * gridsize[d]^-1 * (n_d-1));
                        G_d   = Vec{N}(recipcell[:, d] * centered_index(n_d-1, gridsize[d]/2));
                        G2_d  = norm(G_d)^2))
    innerloop = :(@nexprs $N d -> ((@nref $N r n)  += vec_d;
                                   (@nref $N G n)  += G_d;
                                   (@nref $N G2 n) += G2_d))
    quote
        @assert size(cell)[1] == N error("Cell and grid dimension need to be the same")
        cell = ustrip.(cell .|> u"a₀")
        recipcell = 2pi * (cell^-1)'
        r  = zeros(Vec{N, T}, gridsize...)
        G  = zeros(Vec{N, T}, gridsize...)
        G2 = zeros(T, gridsize...)
        Ω  = det(cell)
        @nloops $N n $dimrange $preloop $innerloop
        return PWGrid(r, G, G2, Ω)
    end
end
Base.size(x::PWGrid) = size(x.r)
Base.length(x::PWGrid) = length(x.r)
Ω(x::PWGrid) = x.Ω
dΩ(x::PWGrid) = Ω(x)/length(x)

Δ!(x, grid::PWGrid) = (x .*= -Ω(grid) .* grid.G2)
Δ(x, grid::PWGrid)  = (t = copy(x); Δ!(t, grid); t)

invΔ!(x, grid::PWGrid{T}) where T = (x[1] = zero(T); x[2:end] ./= (-Ω(grid) .* grid.G2[2:end]))
invΔ(x, grid::PWGrid)  = (t = copy(x); invΔ!(t, grid); t)

function structure_factor(crystal, grid::PWGrid{T}) where T
    G  = grid.G
    Sf = zeros(Complex{T}, size(G))

    tdot = zero(T)
    for atom in eachatom(crystal)
        pos = ustrip.(atom[:position] .|> a₀)
        Sf .+= atom[:charge] .* exp.(-1im .* dot.((pos,), G))
    end
    return Sf
end

function solve_poisson(ρ, grid::PWGrid)
    ρk = -4π .* fft(ρ)
    invΔ!(ρk, grid)
    real.( ifft!(ρk))
end
