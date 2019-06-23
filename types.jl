import Unitful: 𝐋, FreeUnits, Quantity
const DirectSpaceType{T,A}     = Quantity{T,𝐋,FreeUnits{A,𝐋,nothing}}
const ReciprocalSpaceType{T,A} = Quantity{T,𝐋^-1, FreeUnits{A,𝐋^-1,nothing}}

const Point{N, F} = SVector{N, F}
const Point3{F} = SVector{3, F}
const Vec3{F}   = SVector{3, F}

const Mat3{F}   = SMatrix{3, 3, F, 9}
const Mat4{F}   = SMatrix{4, 4, F, 16}

@inline function StaticArrays._inv(::StaticArrays.Size{(3,3)}, A::SMatrix{3,3, LT}) where {LT<:Length}

    @inbounds x0 = SVector{3}(A[1], A[2], A[3])
    @inbounds x1 = SVector{3}(A[4], A[5], A[6])
    @inbounds x2 = SVector{3}(A[7], A[8], A[9])

    y0 = cross(x1,x2)
    d  = StaticArrays.bilinear_vecdot(x0, y0)
    x0 = x0 / d
    y0 = y0 / d
    y1 = cross(x2,x0)
    y2 = cross(x0,x1)

    @inbounds return SMatrix{3, 3}((y0[1], y1[1], y2[1], y0[2], y1[2], y2[2], y0[3], y1[3], y2[3]))
end

@with_kw struct Cell{T<:AbstractFloat,DT<:DirectSpaceType{T},RT<:ReciprocalSpaceType{T}}
	real ::Mat3{DT}
	recip::Mat3{RT} = 2π * inv(real)'
end

abstract type ComponentData end
abstract type AbstractComponent{T <: ComponentData} end
Base.eltype(::AbstractComponent{T}) where {T <: ComponentData} = T

struct Component{T <: ComponentData} <: AbstractComponent{T}
	id  ::Int
	data::Vector{T}
end

struct Atom
	id ::Int
end

struct DFTStructure{T<:AbstractFloat,CT<:Cell{T}}
	cell ::CT
	atoms::Vector{Int}
end

