export SchobelZhuParams, evaluateCharFunc, evaluateLogCharFunc, computeCumulants

struct SchobelZhuParams{T}
    v0::T
    κ::T
    θ::T
    ρ::T
    σ::T
end

#CC the Nemo arbField or Complex type. z: Nemo complex number or Complex number.
function evaluateCharFunc(CC, p::SchobelZhuParams, z::CT, τ::T) where {T,CT}
    return exp(evaluateLogCharFunc(CC, p, z, τ))
end

function evaluateLogCharFunc(CC, p::SchobelZhuParams, z::CT, τ::T) where {T,CT}
    #Lord and Kahl formulation
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = CC(0, 1)
    α = -(z * (z + cc1)) / 2
    β = 2 * (κ - cc1 * σ * ρ * z)
    γ = 2 * σ^2
    D = sqrt(β^2 - 4 * α * γ)
    G = (β - D) / (β + D)
    eDT = exp(-D * τ)
    l = log((1 - G * eDT) / (1 - G))
    sqrteDT = sqrt(eDT)
    As =
        (β - D) * κ^2 * θ^2 / (2 * D^3 * σ^2) * (
            β * (D * τ - 4) +
            D * (D * τ - 2) +
            4 * sqrteDT * ((D^2 - 2 * β^2) / (β + D) * sqrteDT + 2 * β) / (1 - G * eDT)
        )
    A = As + (β - D) * τ / 4 - l / 2
    Bs = 2 * κ * θ * (β - D) * (1 - sqrteDT)^2 / (D * γ * (1 - G * eDT))
    Bv = (1 - eDT) / (1 - G * eDT) * (β - D) / (2 * γ)
    return A + Bs * v0 + Bv * v0^2
end

function evaluateLogCharFuncSZ(CC, p::SchobelZhuParams, z::CT, τ::T) where {T,CT}
    #Zhu formulation
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = CC(0, 1)
    iu = complex(0, 1) * z
    s12 = -iu / 2 * (iu * (1 - ρ^2) - 1 + 2 * ρ * κ / σ)
    s22 = ρ * κ * θ / σ * iu
    s32 = ρ / (2 * σ) * iu
    γ1 = sqrt(2 * σ^2 * s12 + κ^2)
    γ3 = κ^2 * θ - s22 * σ^2
    γ2 = (κ - 2 * σ^2 * s32) / γ1
    ch1 = cosh(γ1 * τ)
    sh1 = sinh(γ1 * τ)
    γ4 = ch1 + γ2 * sh1
    H3 = (κ - γ1 * (sh1 + γ2 * ch1) / γ4) / σ^2
    H4 = ((κ * θ * γ1 - γ2 * γ3) + γ3 * (sh1 + γ2 * ch1)) / (γ4 * γ1 * σ^2) - κ * θ / σ^2
    H5 =
        -log(γ4) / 2 +
        sh1 / (2 * γ1^3 * γ4 * σ^2) * ((κ * θ * γ1 - γ2 * γ3)^2 - γ3^2 * (1 - γ2^2))
    H5 +=
        γ3 / (γ1^3 * γ4 * σ^2) * (κ * θ * γ1 - γ2 * γ3) * (γ4 - 1) +
        τ / (2 * γ1^2 * σ^2) * (κ * γ1^2 * (σ^2 - κ * θ^2) + γ3^2)
    return -s32 * v0^2 - iu * ρ * σ * τ / 2 + H3 / 2 * v0^2 + H4 * v0 + H5
end

function computeCumulants(p::SchobelZhuParams{T}, τ::T) where {T}
    κ = p.κ
    θ = p.θ
    v0 = p.v0
    σ = p.σ
    ρ = p.ρ

    shk = sinh(κ * τ)
    chk = cosh(κ * τ)
    c1H3 = (chk * (chk - shk) * σ + ρ * κ - σ) / (σ * κ)
    c1H4 = (chk - 1) * θ * (shk - chk) / κ
    c1H5 =
        (
            (2 * κ * θ^2 - σ^2) * chk^2 +
            ((-2 * κ * θ^2 + σ^2) * shk - 4 * κ * θ^2) * chk +
            4 * shk * κ * θ^2 +
            2 * τ * (ρ * σ - θ^2) * κ^2 +
            (-σ^2 * τ + 2 * θ^2) * κ +
            σ^2
        ) / (4 * κ^2)
    c1 = -ρ / (2 * σ) * v0^2 - ρ * σ * τ / 2 + c1H3 / 2 * v0^2 + c1H4 * v0 + c1H5
    c2H3 =
        (
            (-2 * κ * ρ * σ + 2 * κ^2) * chk^2 -
            2 * (κ * ρ * σ - 1 / 2 * σ^2 - κ^2) * shk * chk +
            (2 * ρ * σ * τ - 2) * κ^2 +
            (-σ^2 * τ + 2 * ρ * σ) * κ
        ) / (κ^3 * (2 * chk * shk + 2 * chk^2 - 1))
    c2H4 =
        2 *
        (
            (-2 * κ * ρ * σ + κ^2 + σ^2) * chk^2 +
            (
                (-2 * κ * ρ * σ + 1 / 2 * σ^2 + κ^2) * shk +
                (ρ * σ * τ - 1) * κ^2 +
                (-1 / 2 * σ^2 * τ + 2 * ρ * σ) * κ - σ^2
            ) * chk +
            ((ρ * σ * τ - 1) * κ^2 + (-1 / 2 * σ^2 * τ + 2 * ρ * σ) * κ - σ^2 / 2) * shk -
            κ * σ * (ρ * κ - σ / 2) * τ
        ) *
        θ / (κ^3 * (2 * chk * shk + 2 * chk^2 - 1))
    c2H5 =
        (
            (-8 * κ * σ^2 * θ^2 + 2 * σ^4) * chk^4 -
            2 * σ^2 * ((-4 * κ * θ^2 + σ^2) * shk - 8 * κ * θ^2) * chk^3 -
            8 *
            κ *
            (
                2 * σ^2 * θ^2 * shk +
                (-2 * ρ * σ * τ * θ^2 + θ^2) * κ^2 +
                σ * (ρ * τ * σ^2 + (τ * θ^2 - 1 / 2) * σ - 3 * ρ * θ^2) * κ +
                σ^2 * (-1 / 2 * σ^2 * τ + ρ * σ + 2 * θ^2)
            ) *
            chk^2 +
            (
                (
                    (-16 * ρ * σ * τ * θ^2 + 8 * θ^2) * κ^3 +
                    (8 * ρ * τ * σ^3 + (8 * τ * θ^2 - 4) * σ^2 - 24 * ρ * σ * θ^2) * κ^2 +
                    (-4 * σ^4 * τ + 8 * ρ * σ^3 + 20 * σ^2 * θ^2) * κ - σ^4
                ) * shk -
                16 *
                ((ρ * σ * τ - 1) * κ^2 + (-1 / 2 * σ^2 * τ + 3 * ρ * σ) * κ - σ^2) *
                κ *
                θ^2
            ) * chk +
            16 *
            κ *
            θ^2 *
            ((ρ * σ * τ - 1) * κ^2 + (-1 / 2 * σ^2 * τ + 3 * ρ * σ) * κ - (3 * σ^2) / 2) *
            shk +
            8 * κ^4 * τ * θ^2 +
            (-24 * ρ * σ * τ * θ^2 + 4 * σ^2 * τ - 8 * θ^2) * κ^3 +
            ((12 * τ * θ^2 - 4) * σ^2 + 24 * ρ * σ * θ^2) * κ^2 +
            (-σ^4 * τ + 8 * ρ * σ^3 - 8 * σ^2 * θ^2) * κ - 2 * σ^4
        ) / (8 * κ^4)
    c2 = c2H3 / 2 * v0^2 + c2H4 * v0 + c2H5

    c4H3 = 32*σ^2*(((3*ρ^2 + 3/4)*κ^4 - 3*σ*ρ*(ρ^2 + 3/2)*κ^3 + 9*(ρ^2 + 1/4)*σ^2*κ^2/2 - (15*κ*ρ*σ^3)/8 + (3*σ^4)/16)*chk^4 - 3*((-ρ^2 - 1/4)*κ^4 + σ*ρ*(ρ^2 + 3/2)*κ^3 - 3*(ρ^2 + 1/4)*σ^2*κ^2/2 + (5*κ*ρ*σ^3)/8 - (3*σ^4)/32)*shk*chk^3 + (ρ^2*τ^2*(-3/2 + ρ*σ*τ)*κ^6 + (3*τ*(-ρ^2*σ^2*τ^2 + ρ*τ*(ρ^2 + 2)*σ - ρ^2 - 1/2)*κ^5)/2 + ((3*ρ*σ^3*τ^3)/4 - 3*τ^2*(ρ^2 + 3/8)*σ^2 + (3*ρ*τ*(ρ^2 + 2)*σ)/2 - 3/4 - (15*ρ^2)/4)*κ^4 + (15*σ*(-σ^3*τ^3/30 + (2*ρ*τ^2*σ^2)/5 + (-3/5*ρ^2*τ - 3/20*τ)*σ + ρ^3 + (6*ρ)/5)*κ^3)/4 - (9*σ^2*(1/24*σ^2*τ^2 - 5/24*ρ*σ*τ + ρ^2 + 1/4)*κ^2)/2 + 15*(-σ*τ/10 + ρ)*σ^3*κ/8 - (3*σ^4)/16)*chk^2 + (ρ^2*τ^2*(-3/2 + ρ*σ*τ)*κ^6 + (3*τ*(-ρ^2*σ^2*τ^2 + ρ*τ*(ρ^2 + 2)*σ - ρ^2 - 1/2)*κ^5)/2 + ((3*ρ*σ^3*τ^3)/4 - 3*τ^2*(ρ^2 + 3/8)*σ^2 + (3*ρ*τ*(ρ^2 + 2)*σ)/2 - 3/8 - (9*ρ^2)/4)*κ^4 + (9*σ*(-σ^3*τ^3/18 + (2*ρ*τ^2*σ^2)/3 + (-ρ^2*τ - 1/4*τ)*σ + ρ^3 + ρ)*κ^3)/4 - 9*(1/12*σ^2*τ^2 - 1/6*ρ*σ*τ + ρ^2 + 1/8)*σ^2*κ^2/4 + (15*σ^3*(σ*τ/5 + ρ)*κ)/32 - (3*σ^4)/64)*shk*chk - ((ρ^2*τ^2*(-3/2 + ρ*σ*τ)*κ^5 + (3*τ*(-ρ^2*σ^2*τ^2 + ρ*τ*(ρ^2 + 2)*σ - ρ^2 - 1/2)*κ^4)/2 + ((3*ρ*σ^3*τ^3)/4 - (3*τ^2*(ρ^2 + 3/4)*σ^2)/2 + (3*ρ*τ*(ρ^2 + 1)*σ)/2 - (3*ρ^2)/2)*κ^3 + (-1/8*σ^4*τ^3 + 3/2*ρ^3*σ + 9/16*σ^2*τ)*κ^2 - (15*σ^3*(-σ*τ/5 + ρ)*τ*κ)/16 + (3*σ^4*τ)/32)*κ)/2)/(κ^7*(8*chk^4 + 8*chk^3*shk - 8*chk^2 - 4*chk*shk + 1))
    c4H4 = -32*σ^2*θ*(-(3*chk^8*σ^4)/4 + (3*σ^4*(shk + 1)*chk^7)/4 + (9*σ^2*(-σ^2*shk/6 + κ*((ρ*σ*τ - 1/2)*κ + σ*(-σ*τ/2 + ρ)))*chk^6)/2 - 9*(((ρ*σ*τ - 1/2)*κ^2 + σ*(-σ*τ/2 + ρ)*κ - σ^2/12)*shk + 5*((ρ*σ*τ - 3/5)*κ - σ^2*τ/2 + (6*ρ*σ)/5)*κ/6)*σ^2*chk^5/2 + (15*((ρ*σ*τ - 3/5)*κ^2 + (-1/2*σ^2*τ + 6/5*ρ*σ)*κ - σ^2/10)*σ^2*shk/4 - 6*κ*((ρ^2*σ^2*τ^2 - ρ*σ*τ + 1/8)*κ^3 + 2*σ*(-ρ*τ^2*σ^2/2 + τ*(ρ^2 + 3/8)*σ - ρ/2)*κ^2 + (1/4*τ^2*σ^4 - 7/8*ρ*τ*σ^3 + 5/4*ρ^2*σ^2)*κ - σ^4*τ/8))*chk^4 + (((6*ρ^2*σ^2*τ^2 - 6*ρ*σ*τ + 3/4)*κ^4 + 12*σ*(-ρ*τ^2*σ^2/2 + τ*(ρ^2 + 3/8)*σ - ρ/2)*κ^3 + (15*σ^2*(1/5*σ^2*τ^2 - ρ*σ*τ + ρ^2 + 3/20)*κ^2)/2 - (9*σ^3*(-σ*τ/6 + ρ)*κ)/4 + (9*σ^4)/32)*shk + (27*κ*((ρ^2*σ^2*τ^2 - 4/3*ρ*σ*τ + 2/9)*κ^3 + (8*σ*(-(3*ρ*τ^2*σ^2)/8 + τ*(ρ^2 + 3/8)*σ - (2*ρ)/3)*κ^2)/3 + (1/4*τ^2*σ^4 - 10/9*ρ*τ*σ^3 + 20/9*ρ^2*σ^2)*κ - (7*σ^4*τ)/36))/8)*chk^3 + (((-27/8*ρ^2*σ^2*τ^2 + 9/2*ρ*σ*τ - 3/4)*κ^4 - 9*σ*(-(3*ρ*τ^2*σ^2)/8 + τ*(ρ^2 + 3/8)*σ - (2*ρ)/3)*κ^3 - (15*σ^2*(9/80*σ^2*τ^2 - 3/4*ρ*σ*τ + ρ^2 + 3/20)*κ^2)/2 + 9*(-σ*τ/8 + ρ)*σ^3*κ/4 - (9*σ^4)/32)*shk + κ*(ρ^2*τ^2*(-3/2 + ρ*σ*τ)*κ^5 + 3*τ*(-ρ^2*σ^2*τ^2/2 + ρ*τ*(ρ^2 + 1)*σ - ρ^2 - 1/4)*κ^4 + ((3*ρ*σ^3*τ^3)/4 - (9*σ^2*τ^2)/8 + 9*(ρ^2 + 1/3)*τ*ρ*σ/2 - (9*ρ^2)/4)*κ^3 + 3*σ*(-σ^3*τ^3/24 - (7*ρ*τ^2*σ^2)/8 + τ*(ρ^2/2 + 9/16)*σ + ρ^3)*κ^2 - 27*(-(5*σ*τ)/18 + ρ)*σ^3*τ*κ/8 + (3*σ^4*τ)/8))*chk^2 + (((-ρ^3*σ*τ^3 + 3/2*ρ^2*τ^2)*κ^6 - 3*τ*(-ρ^2*σ^2*τ^2/2 + ρ*τ*(ρ^2 + 1)*σ - ρ^2 - 1/4)*κ^5 + (-(3*ρ*σ^3*τ^3)/4 + 3*τ^2*(ρ^2 + 3/8)*σ^2 - (9*ρ*τ*(ρ^2 + 1)*σ)/2 + 3/8 + (9*ρ^2)/4)*κ^4 - 3*(-σ^3*τ^3/24 + ρ*τ^2*σ^2/8 + τ*(-(3*ρ^2)/2 - 3/16)*σ + ρ^3 + ρ)*σ*κ^3 + 15*(-1/20*σ^2*τ^2 - 1/4*ρ*σ*τ + ρ^2 + 9/40)*σ^2*κ^2/4 - (27*σ^3*(-σ*τ/18 + ρ)*κ)/16 + (15*σ^4)/64)*shk - κ*((ρ^3*σ*τ^3 - 3*ρ^2*τ^2)*κ^5 + 6*(-ρ^2*σ^2*τ^2/4 + ρ*τ*(ρ^2 + 1)*σ - 2*ρ^2 - 1/2)*τ*κ^4 + ((3*ρ*σ^3*τ^3)/4 + τ^2*(9*ρ^2 - 9/4)*σ^2 + 18*ρ^3*σ*τ - 18*ρ^2)*κ^3 + (-σ^4*τ^3/8 - (57*ρ*τ^2*σ^3)/4 + τ*(18*ρ^2 + 45/4)*σ^2 + 24*ρ^3*σ)*κ^2 + (33/8*τ^2*σ^4 - 45/2*ρ*τ*σ^3)*κ + (21*σ^4*τ)/8)/8)*chk + (ρ^2*τ^2*(ρ*σ*τ - 3)*κ^6/8 + 3*(-ρ^2*σ^2*τ^2/4 + ρ*τ*(ρ^2 + 1)*σ - 2*ρ^2 - 1/2)*τ*κ^5/4 + ((3*ρ*σ^3*τ^3)/32 - 9*(ρ^2 + 1/2)*τ^2*σ^2/16 + (9*ρ*τ*(ρ^2 + 1)*σ)/4 - 3/8 - (9*ρ^2)/4)*κ^4 + 3*σ*(-σ^3*τ^3/192 - ρ*τ^2*σ^2/32 + τ*(-(3*ρ^2)/4 - 3/32)*σ + ρ^3 + ρ)*κ^3 - (15*σ^2*(-1/40*σ^2*τ^2 - 1/8*ρ*σ*τ + ρ^2 + 9/40)*κ^2)/4 + 27*(-σ*τ/36 + ρ)*σ^3*κ/16 - (15*σ^4)/64)*shk - κ*(ρ^2*τ*(-3/2 + ρ*σ*τ)*κ^5 + (-(3*ρ^2*σ^2*τ^2)/2 + 3*ρ*τ*(ρ^2 + 1)*σ - 3*ρ^2 - 3/4)*κ^4 + (9*σ*(ρ*τ^2*σ^2/6 + τ*(-ρ^2 - 1/4)*σ + ρ^3 + (4*ρ)/3)*κ^3)/2 - 15*(1/60*σ^2*τ^2 - 1/4*ρ*σ*τ + ρ^2 + 9/40)*σ^2*κ^2/2 + (27*σ^3*(-σ*τ/18 + ρ)*κ)/8 - (15*σ^4)/32)*τ/2)/κ^7
    c4H5=-8*σ^2*(-(3*σ^4*(-8*κ*θ^2 + σ^2)*chk^8)/16 + (3*σ^4*((-8*κ*θ^2 + σ^2)*shk - 16*κ*θ^2)*chk^7)/16 + (3*κ*σ^2*(2*σ^2*θ^2*shk + (-6*ρ*σ*τ*θ^2 + 3*θ^2)*κ^2 + σ*(ρ*τ*σ^2 + (3*τ*θ^2 - 1/2)*σ - 7*ρ*θ^2)*κ + σ^2*(-1/2*σ^2*τ + ρ*σ + 2*θ^2))*chk^6)/2 - (3*σ^2*(((-6*ρ*σ*τ*θ^2 + 3*θ^2)*κ^3 + σ*(ρ*τ*σ^2 + (3*τ*θ^2 - 1/2)*σ - 7*ρ*θ^2)*κ^2 + σ^2*(-1/2*σ^2*τ + ρ*σ + 5/2*θ^2)*κ - σ^4/16)*shk - 10*((ρ*σ*τ - 3/5)*κ^2 + (-1/2*σ^2*τ + 7/5*ρ*σ)*κ - σ^2/5)*κ*θ^2)*chk^5)/2 - 3*κ*(5*σ^2*((ρ*σ*τ - 3/5)*κ^2 + (-1/2*σ^2*τ + 7/5*ρ*σ)*κ - (3*σ^2)/10)*θ^2*shk - 4*(ρ^2*σ^2*τ^2 - ρ*σ*τ + 1/8)*θ^2*κ^4 + σ*(ρ^2*σ^3*τ^2 + (4*ρ*τ^2*θ^2 - ρ*τ)*σ^2 + (1/8 + (-10*ρ^2 - 3)*τ*θ^2)*σ + 5*ρ*θ^2)*κ^3 + 2*(-ρ*τ^2*σ^3/2 + τ*(-τ*θ^2/2 + ρ^2 + 3/8)*σ^2 + (17/4*ρ*τ*θ^2 - 1/2*ρ)*σ - 15*(ρ^2 + 2/5)*θ^2/4)*σ^2*κ^2 + (σ^6*τ^2/4 - ρ*σ^5*τ + (-(3*τ*θ^2)/2 + (5*ρ^2)/4)*σ^4 + 7*ρ*σ^3*θ^2)*κ - θ^2*σ^4 - σ^6*τ/16)*chk^4 + ((-12*(ρ^2*σ^2*τ^2 - ρ*σ*τ + 1/8)*θ^2*κ^5 + 3*σ*(ρ^2*σ^3*τ^2 + (4*ρ*τ^2*θ^2 - ρ*τ)*σ^2 + (1/8 + (-10*ρ^2 - 3)*τ*θ^2)*σ + 5*ρ*θ^2)*κ^4 + 6*σ^2*(-ρ*τ^2*σ^3/2 + τ*(-τ*θ^2/2 + ρ^2 + 3/8)*σ^2 + (5*ρ*τ*θ^2 - 1/2*ρ)*σ - 15*(ρ^2 + 1/2)*θ^2/4)*κ^3 + 15*(σ^3*τ^2/5 - ρ*τ*σ^2 + (-(9*τ*θ^2)/5 + ρ^2 + 1/10)*σ + 7*ρ*θ^2)*σ^3*κ^2/4 - (3*σ^4*(-1/4*σ^2*τ + ρ*σ + 27/4*θ^2)*κ)/4 + (9*σ^6)/128)*shk - (27*κ*((ρ^2*σ^2*τ^2 - 4/3*ρ*σ*τ + 2/9)*κ^4 + (-ρ*τ^2*σ^3 + τ*((10*ρ^2)/3 + 1)*σ^2 - (20*ρ*σ)/9)*κ^3 + (10*σ^2*(3/40*σ^2*τ^2 - 19/30*ρ*σ*τ + ρ^2 + 1/5)*κ^2)/3 - (14*σ^3*(-(5*σ*τ)/56 + ρ)*κ)/9 + (2*σ^4)/9)*θ^2)/2)*chk^3 + κ*((27*θ^2*((ρ^2*σ^2*τ^2 - 4/3*ρ*σ*τ + 2/9)*κ^4 + (-ρ*τ^2*σ^3 + τ*((10*ρ^2)/3 + 1)*σ^2 - (20*ρ*σ)/9)*κ^3 + (τ^2*σ^4/4 - (8*ρ*τ*σ^3)/3 + ((10*ρ^2)/3 + 1)*σ^2)*κ^2 - (7*σ^3*(-(5*σ*τ)/28 + ρ)*κ)/3 + (5*σ^4)/12)*shk)/2 + (-2*ρ^3*σ*τ^3*θ^2 + 3*ρ^2*τ^2*θ^2)*κ^6 + (ρ^3*σ^3*τ^2 + 3*τ*ρ^2*(τ*θ^2 - 1/2)*σ^2 - 9*(ρ^2 + 2/3)*τ*ρ*θ^2*σ + 9*(ρ^2 + 1/6)*θ^2)*τ*κ^5 + (-(3*ρ^2*σ^4*τ^3)/2 + 3*(-τ*θ^2/2 + ρ^2 + 1)*τ^2*ρ*σ^3 + 9*(ρ^2 + 1/4)*τ*(τ*θ^2 - 1/3)*σ^2 - 18*ρ*τ*θ^2*(ρ^2 + 1)*σ + (9*ρ^2 + 3)*θ^2)*κ^4 + (9*σ*(ρ*σ^4*τ^3/6 - ((-τ*θ^2/9 + ρ^2 + 1/2)*τ^2*σ^3)/2 + ρ*τ*(-τ*θ^2/3 + ρ^2 + 2/3)*σ^2 + (τ*((20*ρ^2)/3 + 5/4)*θ^2 - ρ^2/2)*σ - (10*ρ*θ^2*(ρ^2 + 2))/3)*κ^3)/2 + 3*σ^2*(-σ^4*τ^3/24 - ρ*τ^2*σ^3/8 - ((τ*θ^2/4 + ρ^2 - 3/8)*τ*σ^2)/2 + (-27/8*ρ*τ*θ^2 + ρ^3)*σ + (15*ρ^2 + 3)*θ^2)*κ^2 + (-21*ρ*σ^3*θ^2 - 9/8*ρ*σ^5*τ + 3/8*σ^6*τ^2)*κ + 3*θ^2*σ^4 + (3*σ^6*τ)/32)*chk^2 + (((2*ρ^3*σ*τ^3*θ^2 - 3*ρ^2*τ^2*θ^2)*κ^7 - (ρ^3*σ^3*τ^2 + 3*τ*ρ^2*(τ*θ^2 - 1/2)*σ^2 - 9*(ρ^2 + 2/3)*τ*ρ*θ^2*σ + 9*(ρ^2 + 1/6)*θ^2)*τ*κ^6 + ((3*ρ^2*σ^4*τ^3)/2 - 3*(-τ*θ^2/2 + ρ^2 + 1)*τ^2*ρ*σ^3 + (-15*τ^2*(ρ^2 + 3/20)*θ^2 + τ*(3*ρ^2 + 3/4))*σ^2 + 18*(ρ^2 + 4/3)*τ*ρ*θ^2*σ - 9*(ρ^2 + 5/12)*θ^2)*κ^5 - (9*σ*(ρ*σ^4*τ^3/6 - (5*τ^2*(-τ*θ^2/15 + ρ^2 + 3/10)*σ^3)/6 + ρ*τ*(-(5*τ*θ^2)/3 + ρ^2 + 1)*σ^2 + (τ*(10*ρ^2 + 9/4)*θ^2 - ρ^2/2 - 1/24)*σ - 10*(ρ^2 + 5/2)*ρ*θ^2/3)*κ^4)/2 - 3*σ^2*(-σ^4*τ^3/24 + (3*ρ*τ^2*σ^3)/8 - 3*(-τ*θ^2/4 + ρ^2 + 1/8)*τ*σ^2/2 + ρ*(-(35*τ*θ^2)/4 + ρ^2 + 1/2)*σ + 75*(ρ^2 + 27/100)*θ^2/4)*κ^3 + (15*σ^3*(-ρ*τ*σ^2/2 + (-(21*τ*θ^2)/10 + ρ^2 + 3/20)*σ + (189*ρ*θ^2)/10)*κ^2)/8 - (9*σ^4*(-1/6*σ^2*τ + ρ*σ + 65/6*θ^2)*κ)/16 + (15*σ^6)/256)*shk + κ*((ρ^3*σ*τ^3 - 3*ρ^2*τ^2)*κ^6 + 9*(-ρ^2*σ^2*τ^2/6 + ρ*τ*(ρ^2 + 2/3)*σ - 2*ρ^2 - 1/3)*τ*κ^5 + ((3*ρ*σ^3*τ^3)/4 + 3*τ^2*(ρ^2 - 3/4)*σ^2 + (36*ρ^3 + 18*ρ)*τ*σ - 36*ρ^2 - 6)*κ^4 + 60*σ*(-σ^3*τ^3/480 - (7*ρ*τ^2*σ^2)/40 + τ*(-ρ^2/4 + 3/80)*σ + ρ^3 + ρ)*κ^3 - 90*σ^2*(-3/80*σ^2*τ^2 + 1/8*ρ*σ*τ + ρ^2 + 1/5)*κ^2 + (27/8*σ^4*τ + 42*ρ*σ^3)*κ - 6*σ^4)*θ^2/2)*chk - κ*((ρ^3*σ*τ^3 - 3*ρ^2*τ^2)*κ^6 + 9*(-ρ^2*σ^2*τ^2/6 + ρ*τ*(ρ^2 + 2/3)*σ - 2*ρ^2 - 1/3)*τ*κ^5 + ((3*ρ*σ^3*τ^3)/4 - 21*(ρ^2 + 3/14)*τ^2*σ^2/2 + 36*ρ*τ*(ρ^2 + 1)*σ - 36*ρ^2 - 9)*κ^4 + 60*(-σ^3*τ^3/480 + ρ*τ^2*σ^2/20 + τ*(-ρ^2 - 3/16)*σ + ρ^3 + (3*ρ)/2)*σ*κ^3 - 135*σ^2*(-19/90*ρ*σ*τ + ρ^2 + 1/4)*κ^2 + (-33/8*σ^4*τ + 315/4*ρ*σ^3)*κ - (105*σ^4)/8)*θ^2*shk/2 + (-3/2 + ρ*σ*τ)*τ^2*ρ^2*θ^2*κ^7 - ((ρ^3*σ^3*τ^2 + 3*τ*ρ^2*(τ*θ^2 - 1/2)*σ^2 - 9*(ρ^2 + 2/3)*τ*ρ*θ^2*σ + 18*(ρ^2 + 1/4)*θ^2)*τ*κ^6)/2 + ((3*ρ^2*σ^4*τ^3)/4 - 3*(-τ*θ^2/2 + ρ^2 + 1)*τ^2*ρ*σ^3/2 - 9*(τ*(ρ^2 + 1/8)*θ^2 - ρ^2/12 - 1/48)*τ*σ^2 + 15*(ρ^2 + 17/10)*τ*ρ*θ^2*σ + 9*(ρ^2 + 1/6)*θ^2)*κ^5 - (3*σ*(ρ*σ^4*τ^3/4 - 3*(-τ*θ^2/18 + ρ^2 + 1/4)*τ^2*σ^3/2 + ρ*τ*(-(7*τ*θ^2)/2 + ρ^2 + 1)*σ^2 + (τ*((55*ρ^2)/2 + 57/8)*θ^2 - (3*ρ^2)/2 - 1/4)*σ + 10*ρ*θ^2*(ρ^2 + 1))*κ^4)/2 - 3*(-σ^4*τ^3/48 + (5*ρ*τ^2*σ^3)/16 - 5*(-τ*θ^2/2 + ρ^2 + 3/20)*τ*σ^2/8 + ρ*(-(141*τ*θ^2)/16 + ρ^2 + 1)*σ - 15*(ρ^2 + 1/5)*θ^2/2)*σ^2*κ^3 + (15*σ^3*(σ^3*τ^2/40 - (3*ρ*τ*σ^2)/20 + (-(49*τ*θ^2)/40 + ρ^2 + 1/5)*σ - (14*ρ*θ^2)/5)*κ^2)/4 - (3*σ^4*(-5/128*σ^2*τ + ρ*σ - θ^2)*κ)/2 + (3*σ^6)/16)/κ^8
    c4 = c4H3 / 2 * v0^2 + c4H4 * v0 + c4H5
    return c1, c2, c4
end