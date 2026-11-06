export CGMYParams, evaluateCharFunc, evaluateLogCharFunc, computeCumulants
import SpecialFunctions: gamma
struct CGMYParams{T}
    C::T
    G::T
    M::T
    Y::T
    σ::T
    ω::T
    CGMYParams(C::T, G::T, M::T, Y::T) where {T} = new{T}(C, G, M, Y, zero(T), driftCorrection(C, G, M, Y))
    CGMYParams(C::T, G::T, M::T, Y::T, σ::T) where {T} = new{T}(C, G, M, Y, σ, driftCorrection(C, G, M, Y))
end

asymptoticLogOrder(model::CGMYParams{T}) where {T} = model.Y

#truncation of the Fourier integral of Filon: log|φ(u)| ~ -σ²τu²/2 - c|u|^Y + d only for |u| well beyond G and M, with a
#constant d that may be large (and singular at Y = 1): the first u, growing geometrically from 10, where the integrand
#|φ(u - i/2)|/(u² + 1/4) is below tol (|φ| decreases in |u|), at most 1/√tol
function computeTruncation(cf::CharFunc{CGMYParams{T1}}, τ::T, tol::T) where {T,T1}
    u, umax = T(10), 1 / sqrt(tol)
    while u < umax && abs(evaluateCharFunc(cf, u - oneim(cf) / 2, τ)) / (u^2 + one(T) / 4) > tol
        u *= T(1.25)
    end
    return min(u, umax)
end

DefaultCharFunc(params::CGMYParams{Float64}) = DefaultCharFunc{CGMYParams{Float64},Complex{Float64}}(params)

@inline function evaluateCharFunc(p::CharFunc{CGMYParams{T},CR}, z::CT, τ::T)::CR where {T,CR,CT}
    return exp(evaluateLogCharFunc(p, z, τ))
end


function evaluateLogCharFunc(cf::CharFunc{CGMYParams{T},CR}, z::CT, τ::T)::CR where {T,CR,CT}
    p = model(cf)
    σ = p.σ
    cc1 = oneim(cf)
    iz = cc1 * z
    q = p.ω
    return iz * q * τ - (z * σ)^2 * τ / 2 + p.C * τ * gamma(-p.Y) * ((p.M - iz)^p.Y - p.M^p.Y + (p.G + iz)^p.Y - p.G^p.Y)
end

function computeCumulants(cf::CharFunc{CGMYParams{T},CR}, τ::T) where {T,CR}
    return computeCumulants(model(cf), τ)
end

function computeCumulants(p::CGMYParams{T}, τ::T) where {T}
    c1 = p.C * τ * gamma(1 - p.Y) * (p.M^(p.Y - 1) - p.G^(p.Y - 1)) 
    c2 = p.σ^2 * τ + p.C * τ * gamma(2 - p.Y) * (p.M^(p.Y - 2) + p.G^(p.Y - 2))
    c4 = p.C * τ * gamma(4 - p.Y) * (p.M^(p.Y - 4) + p.G^(p.Y - 4))
    return c1, c2, c4
end

driftCorrection(p::CGMYParams{T}) where {T} = driftCorrection(p.C, p.G, p.M, p.Y)
driftCorrection(C::T, G::T, M::T, Y::T) where {T} = -C * gamma(-Y) * ((M - one(T))^Y - M^(Y) + (G + one(T))^Y - G^Y)