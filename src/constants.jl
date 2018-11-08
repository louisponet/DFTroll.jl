
import Unitful
import Unitful: @unit, @u_str, me

@unit e₀ "eₒ" ElementaryCharge 1.602176620898e-19*u"C" false
@unit kₑ "kₑ" CoulombForceConstant 1/(4π)u"ϵ0" false
@unit ħ "ħ" ReducedPlanckConstant Unitful.ħ false
@unit a₀ "a₀" BohrRadius 1ħ^2/(1kₑ*me*e₀^2) false
@unit Eₕ "Eₕ" HartreeEnergy me*e₀^4*kₑ^2/(1ħ^2) true
@unit Ry "Ry" RydbergEnergy 0.5Eₕ true
@unit rₑ "rₑ" ClassicalElectronRadius 1e₀^2*kₑ/(1me*Unitful.c^2) false
const α = 1e₀^2*1kₑ/(1Unitful.c*ħ)
const mₚ = 1836.15me
const μ_b = e₀*ħ/(2me)
const ϵ₀ = 1/(4π*kₑ)

@unit ρ         "ρ"         Density                         (1a₀)^-3       false
@unit σₑ        "σₑ"        ContractedDensityGradient       (1a₀)^-8       false
@unit ∂ϵ_∂ρ     "∂ϵ_∂ρ"     FirstDensityDerivative          (1Eₕ)*(1a₀)^3  false
@unit ∂ϵ_∂σ     "∂ϵ_∂σ"     FirstGradientDerivative         (1Eₕ)*(1a₀)^8  false
@unit ∂²ϵ_∂ρ²   "∂²ϵ_∂ρ²"   SecondDensityDerivative         (1Eₕ)*(1a₀)^6  false
@unit ∂²ϵ_∂σ²   "∂²ϵ_∂σ²"   SecondGradientDerivative        (1Eₕ)*(1a₀)^16 false
@unit ∂²ϵ_∂ρ∂σ  "∂²ϵ_∂ρ∂σ"  SecondDensityGradientDerivative (1Eₕ)*(1a₀)^11 false
@unit ∂³ϵ_∂ρ³   "∂³ϵ_∂ρ³"   ThirdDensityDerivative          (1Eₕ)*(1a₀)^9  false
@unit ∂³ϵ_∂σ³   "∂³ϵ_∂σ³"   ThirdGradientDerivative         (1Eₕ)*(1a₀)^24 false
@unit ∂³ϵ_∂ρ²∂σ "∂³ϵ_∂ρ²∂σ" ThirdDensity2GradientDerivative (1Eₕ)*(1a₀)^14 false
@unit ∂³ϵ_∂ρ∂σ² "∂³ϵ_∂ρ∂σ²" ThirdDensityGradient2Derivative (1Eₕ)*(1a₀)^19 false
const ϵ = Eₕ
const rho = ρ
const sig = σₑ
const Eh = Eₕ


const VecArray{T, DIM} = Array{Vec{DIM, T}, DIM}
