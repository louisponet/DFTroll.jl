#based on Abinit functionals
"""
    exchange_correlation(ρ, corr_flavor=1, corr_flag=1)

Calculates the exchange and correlation functionals.
Exchange : Slater
Correlation : LDA parametrization from Monte Carlo data
    `corr_flavor == :PZ` : Based on Perdew-Zunger interpolation
        `corr_flag == 1` : J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
        `corr_flag == 2` : G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
    `corr_flavor == :PW` : Based on Perdew-Wang interpolation
        `corr_flag == 1` : J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
        `corr_flag == 2` : G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
"""
function exchange_correlation(ρ, corr_flavor=:PZ, corr_flag=1)
    if ρ <= 1e-10
        return (εₓ = 0.0, Vₓ = 0.0, εc = 0.0, Vc = 0.0)
    end
    rs = (3/4π/ρ)^(1/3)
    corr = corr_flavor == :PZ ? pz(rs, corr_flag) : pw(rs, corr_flag)
    return (slater(rs)..., corr...)
end

slater(rs) = (εₓ = -0.687247939924714/rs * 2/3, Vₓ = 8/9 * -0.687247939924714 / rs)

function perdew_zunger(rs, corr_flag)
    if rs < 1.0
        a = corr_flag == 1 ?  0.0311 :  0.031091
        b = corr_flag == 1 ? -0.048  : -0.046644
        c = corr_flag == 1 ?  0.0020 :  0.00419
        d = corr_flag == 1 ? -0.0116 : -0.00983
        lnrs = log(rs)
        εc = a*lnrs + b + c*rs*lnrs + d*rs
        Vc = a*lnrs + (b-a)/3 + 2/3*c*rs*ln⁠rs + (2d-c)/3*rs
    else
        a = corr_flag == 1 ? -0.1423 : -0.103756
        b = corr_flag == 1 ?  1.0529 :  0.56371
        c = corr_flag == 1 ?  0.3334 :  0.27358
        rs12 = √rs
        εc = a / (1 + b*rs12 + c*rs)
        Vc = εc^2 * (1 + 7/6*b*rs12 + 4/3*c*rs) / a
    end
    return (εc = εc, Vc = Vc)
end

function perdew_wang(rs, corr_flag)
    a = 0.031091
    if rs < 1 && corr_flag == 2
        c0 = a
        c1 = 0.046644
        c2 = 0.00664
        c3 = 0.01043
        lnrs = log(rs)
        εc = c0*lnrs - c1 + c2*rs*lnrs - c3*rs
        Vc = c0*lnrs - (c1 + c0/3) + 2/3*c2*rs*lnrs - (2c3 + c2)/3*rs
    elseif rs > 100 && corr_flag == 2
        rs32 = rs^(3/2)
        εc = -0.4335/rs + 1.4408/rs32
        Vc = -4/3 * -0.4335/rs + 1.5*1.4408/rs32
    else
        a1 = corr_flag == 1 ? 0.21370 : 0.026481
        b1 = 7.5957
        b2 = 3.5876
        b3 = corr_flag == 1 ? 1.6382  : -0.46647
        b4 = corr_flag == 1 ? 0.49294 : 0.13354
        rs12 = √rs
        rs32 = rs * rs12
        rs2  = rs^2
        om  = 2a*(b1*rs12 + b2*rs + b3*rs32 + b4*rs2)
        dom = 2a*(0.5b1*rs12 + b2*rs + 1.5b3*rs32 + 2b4*rs2)
        olog = log(1 + 1/om)
        εc = -2a * (1 + a1*rs) * olog
        Vc = -2a*olog*(1 + 2/3*a1*rs) - 2/3*a*(1 + a1*rs)*dom / (om*(om + 1))
    end
    return (εc = εc, Vc = Vc)
end
