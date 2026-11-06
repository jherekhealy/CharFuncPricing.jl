export KouParams, HestonKouParams
#Double exponential jump diffusion of Kou (2002): log-jumps J of density pUp η₁ exp(-η₁ x) for x > 0 and (1-pUp) η₂ exp(η₂ x)
#for x < 0, at the rate λ, compensated so that the forward is a martingale (η₁ > 1). KouParams adds a constant volatility σ,
#HestonKouParams the stochastic variance of Heston, independent of the jumps.
struct KouParams{T}
    σ::T
    λ::T
    pUp::T
    η1::T
    η2::T
    function KouParams(σ::T, λ::T, pUp::T, η1::T, η2::T) where {T}
        η1 > 1 || throw(DomainError(η1, "η1 must be larger than 1 for a finite forward"))
        new{T}(σ, λ, pUp, η1, η2)
    end
end

struct HestonKouParams{T}
    heston::HestonParams{T}
    λ::T
    pUp::T
    η1::T
    η2::T
    function HestonKouParams(heston::HestonParams{T}, λ::T, pUp::T, η1::T, η2::T) where {T}
        η1 > 1 || throw(DomainError(η1, "η1 must be larger than 1 for a finite forward"))
        new{T}(heston, λ, pUp, η1, η2)
    end
end

const KouJumps{T} = Union{KouParams{T},HestonKouParams{T}}

DefaultCharFunc(params::KouParams{Float64}) = DefaultCharFunc{KouParams{Float64},Complex{Float64}}(params)
DefaultCharFunc(params::HestonKouParams{Float64}) = DefaultCharFunc{HestonKouParams{Float64},Complex{Float64}}(params)

#E[exp(iz J)] - 1 - iz (E[exp(J)] - 1): the compensated jump exponent per unit of intensity
@inline function kouJumpExponent(p::KouJumps, iz)
    q = 1 - p.pUp
    ζ = p.pUp * p.η1 / (p.η1 - 1) + q * p.η2 / (p.η2 + 1) - 1
    return p.pUp * p.η1 / (p.η1 - iz) + q * p.η2 / (p.η2 + iz) - 1 - iz * ζ
end

function evaluateLogCharFunc(cf::CharFunc{KouParams{T},CR}, z::CT, τ::T) where {T,CR,CT}
    p = model(cf)
    iz = oneim(cf) * z
    return τ * (p.σ^2 * iz * (iz - 1) / 2 + p.λ * kouJumpExponent(p, iz))
end

function evaluateLogCharFunc(cf::CharFunc{HestonKouParams{T},CR}, z::CT, τ::T) where {T,CR,CT}
    p = model(cf)
    iz = oneim(cf) * z
    return evaluateLogCharFunc(DefaultCharFunc{HestonParams{T},CR}(p.heston), z, τ) + τ * p.λ * kouJumpExponent(p, iz)
end

#cumulants c1, c2, c4 of the compensated jumps over τ, from the moments E[Jⁿ] = n! (pUp/η₁ⁿ + (-1)ⁿ (1-pUp)/η₂ⁿ)
function kouJumpCumulants(p::KouJumps{T}, τ::T) where {T}
    q = 1 - p.pUp
    ζ = p.pUp * p.η1 / (p.η1 - 1) + q * p.η2 / (p.η2 + 1) - 1
    m1 = p.pUp / p.η1 - q / p.η2
    m2 = 2 * (p.pUp / p.η1^2 + q / p.η2^2)
    m4 = 24 * (p.pUp / p.η1^4 + q / p.η2^4)
    return p.λ * τ * (m1 - ζ), p.λ * τ * m2, p.λ * τ * m4
end

computeCumulants(cf::CharFunc{KouParams{T},CR}, τ::T) where {T,CR} = computeCumulants(model(cf), τ)
function computeCumulants(p::KouParams{T}, τ::T) where {T}
    j1, j2, j4 = kouJumpCumulants(p, τ)
    return j1 - p.σ^2 * τ / 2, j2 + p.σ^2 * τ, j4
end

computeCumulants(cf::CharFunc{HestonKouParams{T},CR}, τ::T) where {T,CR} = computeCumulants(model(cf), τ)
function computeCumulants(p::HestonKouParams{T}, τ::T) where {T}
    h1, h2, h4 = computeCumulants(p.heston, τ)
    j1, j2, j4 = kouJumpCumulants(p, τ)
    return h1 + j1, h2 + j2, h4 + j4
end

#truncation of the Fourier integral of Filon: the jump part of |φ| stays bounded, the diffusion makes it decay as exp(-σ²τu²/2)
#(without diffusion, the integrand decays as 1/u²)
function computeTruncation(cf::CharFunc{KouParams{T1}}, τ::T, tol::T) where {T,T1}
    p = model(cf)
    uJump = 1 / sqrt(tol)
    return p.σ > 0 ? min(max(sqrt(-2 * log(tol) / (p.σ^2 * τ)), T(10)), uJump) : uJump
end
computeTruncation(cf::CharFunc{HestonKouParams{T1}}, τ::T, tol::T) where {T,T1} = computeTruncation(model(cf).heston, τ, tol)

asymptoticLogOrder(model::KouParams{T}) where {T} = 2 * one(T) #the diffusion dominates: log|φ(u)| ~ -σ²u²τ/2
asymptoticLogOrder(model::HestonKouParams{T}) where {T} = one(T)
