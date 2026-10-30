import SpecialFunctions: gamma
import FFTW: idct

export RoughHestonParams, PadeCharFunc, AdamsCharFunc, RoughHestonCVCharFunc, FejerQuadrature

struct RoughHestonParams{T}
    H::T
    ρ::T
    ν::T
    varianceCurve # a function of time.
end

#Fejer first rule (Chebyshev nodes of the first kind) on [-1,1].
struct FejerQuadrature{T}
    x::Vector{T}
    w::Vector{T}
end

function FejerQuadrature{T}(N::Int) where {T}
    x = @. (cos(((2:2:2N) - 1) * T(pi) / 2N))
    n = N - 1
    m = vcat(sqrt(T(2)), zero(T), [T((1 + (-1)^k) / (1 - k^2)) for k = 2:n])
    w = sqrt(T(2) / (n + 1)) * idct(m)
    return FejerQuadrature{T}(x, w)
end
FejerQuadrature(N::Int) = FejerQuadrature{Float64}(N)

#q is any quadrature rule on [-1,1] with fields x (abscissae) and w (weights).
@inline function integrateRule(q, f, a::T, b::T) where {T}
    bma2 = (b - a) / 2
    bpa2 = (b + a) / 2
    return bma2 * sum(q.w .* f.(bma2 .* q.x .+ bpa2))
end

struct PadeCharFunc{MT,CR,ORDER,Q} <: CharFunc{MT,CR} #model type, return type (e.g. Complex or acb)
    model::MT
    q::Q
end

PadeCharFunc{MT,CR,ORDER}(model::MT, q::Q) where {MT,CR,ORDER,Q} = PadeCharFunc{MT,CR,ORDER,Q}(model, q)
PadeCharFunc{CR}(model::MT, q=FejerQuadrature{Float64}(24)) where {MT,CR} = PadeCharFunc{MT,CR,3}(model, q) #
PadeCharFunc(model::RoughHestonParams{T}) where {T} = PadeCharFunc{Complex{T}}(model)


@inline model(cf::PadeCharFunc) = cf.model

@inline function evaluateCharFunc(p::CharFunc{RoughHestonParams{T},CR}, z::CT, τ::T) where {T,CR,CT}
    return exp(evaluateLogCharFunc(p, z, τ))
end

function evaluateLogCharFunc(p::PadeCharFunc{RoughHestonParams{T},CR,3}, z::CT, τ::T) where {T,CR,CT}
    #integrate pade33.
    #φ(0) = φ(-i) = 1 (martingale): the Padé coefficients are 0/0 there
    if abs(z * (z + 1im)) < sqrt(eps(T))
        return zero(typeof(z * one(T)))
    end
    params = model(p)
    numal = params.ν^(1 / (params.H + 0.5))
    integrand = function (u)
        x = u * numal
        pade(p, z, x)[2] * params.varianceCurve(τ - u)
    end

    return integrateRule(p.q, integrand, zero(T), τ)
end

function pade(p::PadeCharFunc{RoughHestonParams{T},CR,3}, a, x) where {T,CR}
    params = model(p)
    H = params.H
    rho = params.ρ
    al = H + 0.5

    aa = sqrt(a * (a + 1im) - rho^2 * a^2)
    rm = -(1im) * rho * a - aa
    rp = -(1im) * rho * a + aa

    b1 = -a * (a + 1im) / (2 * gamma(1 + al))
    b2 = (1 - a * 1im) * a^2 * rho / (2 * gamma(1 + 2 * al))
    b3 = gamma(1 + 2 * al) / gamma(1 + 3 * al) * (a^2 * (1im + a)^2 / (8 * gamma(1 + al)^2) + (a + 1im) * a^3 * rho^2 / (2 * gamma(1 + 2 * al)))

    g0 = rm
    g1 = -rm / (aa * gamma(1 - al))
    g2 = rm / aa^2 / gamma(1 - 2 * al) * (1 + rm / (2 * aa) * gamma(1 - 2 * al) / gamma(1 - al)^2)

    den = g0^3 + 2 * b1 * g0 * g1 - b2 * g1^2 + b1^2 * g2 + b2 * g0 * g2

    p1 = b1
    p2 = (b1^2 * g0^2 + b2 * g0^3 + b1^3 * g1 + b1 * b2 * g0 * g1 - b2^2 * g1^2 + b1 * b3 * g1^2 + b2^2 * g0 * g2 - b1 * b3 * g0 * g2) / den
    q1 = (b1 * g0^2 + b1^2 * g1 - b2 * g0 * g1 + b3 * g1^2 - b1 * b2 * g2 - b3 * g0 * g2) / den
    q2 = (b1^2 * g0 + b2 * g0^2 - b1 * b2 * g1 - b3 * g0 * g1 + b2^2 * g2 - b1 * b3 * g2) / den
    q3 = (b1^3 + 2 * b1 * b2 * g0 + b3 * g0^2 - b2^2 * g1 + b1 * b3 * g1) / den
    p3 = g0 * q3

    y = x^al
    pade = (p1 * y + p2 * y^2 + p3 * y^3) / (1 + q1 * y + q2 * y^2 + q3 * y^3)
    return pade, (pade - rm) * (pade - rp) / 2
end



## Fractional Adams scheme of Diethelm, Ford and Freed, based on the R code from Jim Gatheral.
#The weights do not depend on u nor on the maturity (up to the factor Δ^α), they are computed once.
#The error is of order 1 + α in the time step. With richardson=true, the log characteristic function is extrapolated
#from n and 2n steps assuming this order, at the cost of n + 2n steps, and of the stability limit in u of the n steps.
struct AdamsWeights{T}
    n::Int
    a0::Vector{T} #corrector weight of the initial value, k = 1..n-1
    aj::Vector{T} #corrector weights, q = 0..n-1
    bj::Vector{T} #predictor weights, r = 0..n
end
function AdamsWeights(α::T, n::Int) where {T}
    g1, g2 = gamma(α + 1), gamma(α + 2)
    a0 = [(k^(α + 1) - (k - α) * (k + 1)^α) / g2 for k = 1:n-1]
    aj = [((q + 2)^(α + 1) + q^(α + 1) - 2 * (q + 1)^(α + 1)) / g2 for q = 0:n-1]
    bj = [((r + 1)^α - r^α) / g1 for r = 0:n]
    return AdamsWeights{T}(n, a0, aj, bj)
end

struct AdamsCharFunc{MT,CR,T} <: CharFunc{MT,CR} #model type, return type (e.g. Complex or acb)
    model::MT
    n::Int
    weights::AdamsWeights{T}
    richardson::Bool
    weights2::AdamsWeights{T} #weights for 2n steps, used with richardson=true
end

function AdamsCharFunc{MT,CR}(model::MT, n::Int; richardson::Bool=false) where {MT<:RoughHestonParams,CR}
    α = model.H + 0.5
    w = AdamsWeights(α, n)
    return AdamsCharFunc{MT,CR,typeof(α)}(model, n, w, richardson, richardson ? AdamsWeights(α, 2n) : w)
end
AdamsCharFunc{CR}(model::MT, n::Int; kwargs...) where {MT,CR} = AdamsCharFunc{MT,CR}(model, n; kwargs...) #
AdamsCharFunc(model::RoughHestonParams{T}, n::Int=200; kwargs...) where {T} = AdamsCharFunc{Complex{T}}(model, n; kwargs...)

@inline model(cf::AdamsCharFunc) = cf.model

RoughHestonCVCharFunc(delegate::Union{AdamsCharFunc{RoughHestonParams{T},CR},PadeCharFunc{RoughHestonParams{T},CR}}) where {T,CR} =
    CVCharFunc{RoughHestonParams{T},BlackParams{T},CR}(
        delegate,
        DefaultCharFunc{BlackParams{T},CR}(
            BlackParams{T}(sqrt(model(delegate).varianceCurve(zero(T))))
        ),
    )



function evaluateLogCharFunc(p::AdamsCharFunc{RoughHestonParams{T},CR}, z::CT, τ::T) where {T,CR,CT}
    #φ(0) = φ(-i) = 1 (martingale)
    if abs(z * (z + 1im)) < sqrt(eps(T))
        return zero(typeof(z * one(T)))
    end
    if p.richardson
        r = 2^(model(p).H + 1.5) #2^(1+α)
        return (r * adamsLogCharFunc(model(p), p.weights2, z, τ) - adamsLogCharFunc(model(p), p.weights, z, τ)) / (r - 1)
    end
    return adamsLogCharFunc(model(p), p.weights, z, τ)
end

function adamsLogCharFunc(params::RoughHestonParams, w::AdamsWeights, z, τ)
    dhA = adams(params, w, z, τ)
    n = w.n
    #trapezoidal rule (a plain sum would make the scheme first order)
    integral = (dhA[1] * params.varianceCurve(τ) + dhA[n+1] * params.varianceCurve(zero(τ))) / 2
    for k = 1:n-1
        integral += dhA[k+1] * params.varianceCurve(τ * (1 - k / n))
    end
    return integral * τ / n
end

#returns the derivative of h on the time grid, the element type follows z (e.g. ForwardDiff duals).
adams(p::AdamsCharFunc, a, bigT) = adams(model(p), p.weights, a, bigT)

function adams(params::RoughHestonParams, w::AdamsWeights, a, bigT)
    rho = params.ρ
    nu = params.ν
    obj(u, h) = -u / 2 * (u + 1im) + rho * nu * u * h * 1im + (nu * h)^2 / 2
    n = w.n
    alpha = params.H + 0.5
    dAlpha = (bigT / n)^alpha
    E = typeof(obj(a, zero(a)))
    hA = zeros(E, n + 1)
    dhA = zeros(E, n + 1) #dhA[k] = obj(a, hA[k]), also the f of the scheme
    a0, aj, bj = w.a0, w.aj, w.bj
    # hA[1] is zero, boundary condition
    dhA[1] = obj(a, hA[1])
    a11 = dAlpha / gamma(alpha + 2)
    hP = dAlpha * bj[1] * dhA[1]
    hA[2] = a11 * alpha * dhA[1] + a11 * obj(a, hP)
    dhA[2] = obj(a, hA[2])
    @inbounds for k = 1:(n-1)
        hP = zero(E)
        for j = 1:k+1
            hP += dhA[j] * bj[k+2-j]
        end
        hC = dhA[1] * a0[k]
        for j = 2:k+1
            hC += dhA[j] * aj[k+2-j]
        end
        hA[k+2] = dAlpha * hC + a11 * obj(a, dAlpha * hP)
        dhA[k+2] = obj(a, hA[k+2])
    end
    return dhA
end
