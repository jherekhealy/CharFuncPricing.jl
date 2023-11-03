export HestonParams, evaluateCharFunc, evaluateLogCharFunc, computeCumulants
export BlackParams, evaluateCharFuncAndDerivative
export CVCharFunc, HestonCVCharFunc, cinf

struct HestonParams{T}
    v0::T
    κ::T
    θ::T
    ρ::T
    σ::T
end

struct BlackParams{T}
    σ::T
end
DefaultCharFunc(params::HestonParams{Float64}) =
    DefaultCharFunc{HestonParams{Float64},Complex}(params)

function evaluateLogCharFunc(
    cf::CharFunc{BlackParams{T},CR},
    z::CT,
    τ::T,
)::CR where {T,CR,CT}
    cc1 = oneim(cf)
    p = cf.model
    return -p.σ^2 * τ / 2 * z * (z + cc1)
end

function evaluateLogCharFuncAndDerivative(
    cf::CharFunc{BlackParams{T},CR},
    z::CT,
    τ::T,
)::Tuple{CR,CR} where {T,CR,CT}
    cc1 = oneim(cf)
    p = cf.model
    phi = -p.σ^2 * τ / 2 * z * (z + cc1)
    phi_d = -p.σ^2 * τ / 2 * (2 * z + cc1)
    return phi, phi_d
end

struct CVCharFunc{MAINT,CONTROLT,CR} <: CharFunc{MAINT,CR}
    main::CharFunc{MAINT,CR}
    control::CharFunc{CONTROLT,CR}
end

model(cf::CVCharFunc) = model(cf.main)
oneim(cf::CVCharFunc) = oneim(cf.main)

HestonCVCharFunc(heston::CharFunc{HestonParams{T},CR}) where {T,CR} =
    CVCharFunc{HestonParams{T},BlackParams{T},CR}(
        heston,
        DefaultCharFunc{BlackParams{T},CR}(
            BlackParams{T}(sqrt(model(heston).v0))
        ),
    )

function evaluateCharFunc(
    p::CVCharFunc{MAINT,CONTROLT,CR},
    z::CT,
    τ::T)::CR where {T,CR,CT,MAINT,CONTROLT}
    phi = evaluateCharFunc(p.main, z, τ)
    phiB = evaluateCharFunc(p.control, z, τ)
    return phi - phiB
end

function evaluateCharFuncAndDerivative(
    p::CVCharFunc{MAINT,CONTROLT,CR},
    z::CT,
    τ::T,
)::Tuple{CR,CR} where {T,CR,CT,MAINT,CONTROLT}
    phi, phi_d = evaluateCharFuncAndDerivative(p.main, z, τ)
    phiB, phiB_d = evaluateCharFuncAndDerivative(p.control, z, τ)
    return (phi - phiB, phi_d - phiB_d)
end

#cinf(params::HestonCVParams{T}, τ::T) where {T} = cinf(params.heston, τ)

cinf(params::HestonParams{T}, τ::T) where {T} =
    (params.v0 + params.κ * params.θ * τ) / params.σ * sqrt(1 - params.ρ^2)

#CC the Nemo arbField or Complex type. z: Nemo complex number or Complex number.
@inline function evaluateCharFunc(
    p::CharFunc{MT,CR},
    z::CT,
    τ::T,
)::CR where {T,CR,CT,MT}
    return exp(evaluateLogCharFunc(p, z, τ))
end

@inline evaluateLogCharFunc(
    p::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CT} = evaluateLogCharFuncAL(p, z, τ)


@inline function evaluateLogCharFuncAL(
    cf::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
)::CR where {T,CR,CT}
    #follows Andersen-Lake fast implementation
    p = model(cf)
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    α = -(z * (z + cc1)) * σ^2
    β = κ - σ * ρ * z * cc1
    D = sqrt(β^2 - α)
    rr = (real(β) * real(D) + imag(β) * imag(D)) > 0 ? α / (β + D) : (β - D)
    y = D == 0 ? -τ / 2 : expm1(-D * τ) / (2 * D)
    l = log1p(-rr * y)
    A = (rr * τ - 2 * l) * (κ * θ / σ^2)
    B = z * (z + cc1) * y / (1 - rr * y)
    return A + B * v0
end

#Faster evaluation of the characteristic function phi(ix) where x is real. Useful for optimal alpha.
@inline function evaluateLogCharFuncAtImag(
    cf::CharFunc{HestonParams{T},CR},
    imz::T,
    τ::T,
)::T where {T,CR}
    #follows Andersen-Lake fast implementation
    p = model(cf)
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    α = (imz * (imz + one(T))) * σ^2
    β = κ + σ * ρ * imz
    D2 = β^2 - α
    local ch, sh
    if D2 >= 0
        D = sqrt(D2)
        ch = cosh(D * τ / 2)
        sh = (D == 0) ? τ / 2 : sinh(D * τ / 2) / D
    else
        D = sqrt(-D2)
        sh, ch = sincos(D * τ / 2)
        sh /= D
    end
    A = κ * θ / σ^2 * (β * τ - log((ch + β * sh)^2))
    B = imz * (imz + 1) * sh / (ch + β * sh)
    return A + B * v0
end

@inline function evaluateCharFuncAndDerivative(
    p::CharFunc{MAINT,CR},
    z::CT,
    τ::T,
)::Tuple{CR,CR} where {T,CR,CT,MAINT}
    arg, arg_d = evaluateLogCharFuncAndDerivative(p, z, τ)
    phi = exp(arg)
    phi_d = arg_d * phi
    return phi, phi_d
end

@inline function evaluateLogCharFuncAndDerivative(
    cf::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
)::Tuple{CR,CR} where {T,CR,CT}
    #derivative towards real part of z
    p = model(cf)
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    α = -(z * (z + cc1)) * σ^2
    α_d = -(2 * z + cc1) * σ^2
    β_d = -σ * ρ * cc1
    β = κ + β_d * z

    D = sqrt(β^2 - α)
    D_d = (2 * β * β_d - α_d) / (2 * D)
    local rr
    local rr_d
    if (real(β) * real(D) + imag(β) * imag(D)) > 0
        rr = α / (β + D)
        rr_d = α_d / (β + D) - rr * (β_d + D_d) / (β + D)
    else
        rr = (β - D)
        rr_d = β_d - D_d
    end
    local y
    local y_d
    if D != 0
        em = expm1(-D * τ)
        y = em / (2 * D)
        y_d = (y * (-2 * D_d * (D * τ + 1)) + (-D_d * τ)) / (2 * D) # ((em+1)*(-D_d * τ)*(2*D) - em*(2*D_d))/(2*D)^2
    else
        y = -τ / 2
        y_d = zero(T)
    end
    l = log1p(-rr * y)
    rry_d = (rr_d * y + rr * y_d)
    l_d = -rry_d / (1 - rr * y)
    A = (rr * τ - 2 * l) * (κ * θ / σ^2)
    A_d = (rr_d * τ - 2 * l_d) * (κ * θ / σ^2)
    Bnum = z * (z + cc1) * y
    Bdenom = (1 - rr * y)
    B = Bnum / Bdenom
    B_d = ((2 * z + cc1) * y + y_d * z * (z + cc1)) / Bdenom + rry_d * B / Bdenom
    return A + B * v0, A_d + B_d * v0
end

struct CuiCharFunc{MT,CR} <: CharFunc{MT,CR} #model type, return type (e.g. Complex or acb),
    model::MT
end

CuiCharFunc(params::HestonParams{Float64}) =
    CuiCharFunc{HestonParams{Float64},Complex,Type}(params, Complex)
@inline model(cf::CuiCharFunc) = cf.model
@inline field(cf::CuiCharFunc) = cf.field
@inline evaluateLogCharFunc(
    p::CuiCharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CF,CT} = evaluateLogCharFuncAL(p, z, τ)


function evaluateLogCharFuncCui(
    cf::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
)::CR where {T,CR,CT}
    p = model(cf)

    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    iu = cc1 * z
    ξ = κ - σ * ρ * iu
    α = iu * (-iu + 1)
    d = sqrt(ξ^2 + σ^2 * α)
    dht = d * (τ / 2)
    sh1 = sinh(dht)
    ch1 = cosh(dht)
    A2v = (d * ch1 + ξ * sh1)
    A1 = α * sh1
    # edt = ch1 - sh1
    # logB = (κ * τ / 2 - dht) - log1p( (ξ - d) * sh1/ d * edt )

    edt = (ch1 + sh1)
    logB = (κ * τ / 2 - dht) - log(A2v / (d * edt))
    return -κ * θ * ρ * τ / σ * iu - A1 / A2v * v0 + 2 * κ * θ / σ^2 * logB
end

function computeCumulants(p::HestonParams{T}, τ::T) where {T}
    lambda = p.κ
    ubar = p.θ
    u0 = p.v0
    eta = p.σ
    rho = p.ρ
    term = τ
    if lambda == 0 && eta == 0
        c1 = -u0 * term / 2
        c2 = u0 * term
        c4 = 0
        return c1, c2, c4
    elseif lambda == 0
        c1 = -u0 * term / 2
        c2 = u0 * term * (1 + term * term * eta * (eta * term / 12 - rho / 2))
        c4 =
            u0 *
            eta^2 *
            term^3 *
            (
                2 * rho^2 - 2 * eta * term * rho +
                eta^4 * term^4 * 17 / 1680 +
                eta^4 * term^5 * rho^2 * 11 / 20 - eta * term * rho^3 / 2 +
                eta^2 * term^2 * 3 / 10 - eta^3 * term^3 * rho * 17 / 120 + 1
            )
        return c1, c2, c4
    end
    elt = exp(-lambda * term)
    c1 = (1 - elt) * (ubar - u0) / (2 * lambda) - ubar * term / 2
    lambda2 = lambda * lambda
    lambda3 = lambda * lambda2
    eta2 = eta * eta
    elt2 = elt * elt
    # correct c2
    c2 =
        1 / (8 * lambda3) * (
            2 *
            u0 *
            (
                lambda2 * 4 * (elt * term * rho * eta + 1 - elt) +
                lambda * (4 * eta * rho * (elt - 1) - 2 * elt * term * eta2) +
                eta2 * (1 - elt2)
            ) +
            ubar * (
                8 * lambda3 * term +
                lambda2 * 8 * (-eta * rho * term * (1 + elt) + elt - 1) +
                lambda * (2 * eta2 * term * (1 + 2 * elt) + eta * rho * 16 * (1 - elt)) +
                eta2 * (-5 + 4 * elt + elt2)
            )
        )
    # paper c2
    # c2 = 1 / (8 * lambda3) * (eta*term*lambda*elt*(u0-ubar)*(8*lambda*rho-4*eta) + lambda*rho*eta*(1-elt)*(16*ubar-8*u0) + 2*ubar*lambda*term*(-4*lambda*rho*eta+eta2+4*lambda2) + eta2*((ubar-2*u0)*elt2+ubar*(6*elt-7)+2*u0) + 8*lambda2*(u0-ubar)*(1-elt))

    lambda4 = lambda2 * lambda2
    lambda5 = lambda4 * lambda
    lambda6 = lambda4 * lambda2
    lambda7 = lambda3 * lambda4
    eta3 = eta * eta2
    term2 = term * term
    term3 = term2 * term
    rho2 = rho * rho
    eta4 = eta2 * eta2
    elt3 = elt2 * elt
    elt4 = elt2 * elt2
    rho3 = rho2 * rho
    c4 =
        -1 / lambda7 *
        2 *
        ubar *
        eta2 *
        (
            (
                (term3 * rho3 * eta - 3 * term2 * rho2) * lambda6 -
                3 *
                term *
                (term2 * rho2 * eta2 - 4 * term * rho * (rho2 + 1) * eta + 8 * rho2 + 2) *
                lambda5 / 2 +
                (
                    3 * term3 * rho * eta3 / 4 - 21 * term2 * (rho2 + 3 / 14) * eta2 / 2 +
                    (18 * term * rho3 + 24 * term * rho) * eta - 18 * rho2 - 3
                ) * lambda4 -
                (
                    term3 * eta3 - 42 * term2 * rho * eta2 +
                    (240 * term * rho2 + 54 * term) * eta - 192 * rho3 - 192 * rho
                ) *
                eta *
                lambda3 / 8 -
                3 *
                eta2 *
                (term2 * eta2 - 35 * term / 2 * rho * eta + 40 * rho2 + 15 / 2) *
                lambda2 / 4 - 27 * eta3 * (term * eta - 20 * rho / 3) * lambda / 16 -
                21 * eta4 / 16
            ) * elt +
            (
                (-3 / 4 + 3 * term * rho * eta - 3 * term2 * rho2 * eta2 / 2) * lambda4 +
                3 *
                eta *
                (term2 * rho * eta2 + (-4 * term * rho2 - 3 * term / 2) * eta + 4 * rho) *
                lambda3 / 2 -
                3 *
                eta2 *
                (term2 * eta2 - 14 * term * rho * eta + 20 * rho2 + 6) *
                lambda2 / 8 + (-15 / 16 * term * eta4 + 9 * eta3 * rho / 2) * lambda -
                21 * eta4 / 32
            ) * elt2 +
            3 *
            eta2 *
            (
                (term * rho * eta - 1) * lambda2 +
                (-term / 2 * eta2 + 2 * rho * eta) * lambda - eta2 / 2
            ) *
            elt3 / 8 - 3 * eta4 * elt4 / 128 +
            (-6 * term * rho2 - 3 * term / 2) * lambda5 +
            ((6 * term * rho3 + 9 * term * rho) * eta + 18 * rho2 + 15 / 4) * lambda4 -
            9 * eta * ((rho2 + 0.25) * term * eta + 8 * rho3 / 3 + 10 * rho / 3) * lambda3 +
            15 * eta2 * (term * rho * eta + 10 * rho2 + 11 / 5) * lambda2 / 4 +
            (-33 / 2 * eta3 * rho - 15 / 32 * term * eta4) * lambda +
            279 * eta4 / 128
        )
    c4 +=
        u0 / lambda7 * (
            2 *
            eta2 *
            (
                (
                    (term3 * rho3 * eta - 3 * term2 * rho2) * lambda6 -
                    3 *
                    term *
                    (
                        term2 * rho2 * eta2 - 2 * term * rho * (rho2 + 2) * eta +
                        4 * rho2 +
                        2
                    ) *
                    lambda5 / 2 +
                    (
                        3 * term3 * rho * eta3 / 4 - 6 * (rho2 + 3 / 8) * term2 * eta2 +
                        6 * term * rho * (rho2 + 2) * eta - 6 * rho2
                    ) * lambda4 -
                    eta *
                    (
                        term3 * eta3 - 24 * term2 * rho * eta2 +
                        (72 * term * rho2 + 18 * term) * eta - 48 * rho3
                    ) *
                    lambda3 / 8 -
                    3 * eta2 * (term2 * eta2 - 7 * term * rho * eta - 3) * lambda2 / 8 -
                    3 * eta3 * (term * eta + 10 * rho) * lambda / 16 + 3 * eta4 / 8
                ) * elt +
                (
                    (6 * term * rho * eta - 3 * term2 * rho2 * eta2 - 3 / 2) * lambda4 +
                    3 *
                    (
                        term2 * rho * eta2 +
                        (-3 * term * rho2 - 3 * term / 2) * eta +
                        3 * rho
                    ) *
                    eta *
                    lambda3 -
                    3 *
                    eta2 *
                    (term2 * eta2 - 10 * term * rho * eta + 12 * rho2 + 3) *
                    lambda2 / 4 - 9 * eta3 * (term * eta - 10 * rho / 3) * lambda / 8 -
                    3 * eta4 / 8
                ) * elt2 +
                9 *
                eta2 *
                (
                    (term * rho * eta - 1) * lambda2 +
                    (-term / 2 * eta2 + 5 / 3 * rho * eta) * lambda - eta2 / 3
                ) *
                elt3 / 8 - 3 * eta4 * elt4 / 32 -
                6 *
                ((rho2 + 1 / 4) * lambda2 - 5 * lambda * rho * eta / 4 + 5 * eta2 / 16) *
                (lambda * rho * eta - eta2 / 4 - lambda2)
            )
        )
    return c1, c2, c4
end
