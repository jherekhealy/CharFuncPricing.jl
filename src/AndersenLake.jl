import Roots: find_zero, A42
import Optim: optimize, Brent, minimizer
import DoubleExponentialFormulas: quadde, quaddeo

export ALCharFuncPricer, TanhSinhQuadrature, ModlobQuadrature, DEQuadrature

abstract type ALQuadrature{T} end


tolerance(q::ALQuadrature) = q.tol
struct ModlobQuadrature{T} <: ALQuadrature{T}
    tol::T
end

struct DEQuadrature{T} <: ALQuadrature{T}
    tol::T
end

struct TanhSinhQuadrature{T} <: ALQuadrature{T}
    h::T
    y::Array{T,1}
    w::Array{T,1}
    tol::T
    function TanhSinhQuadrature(n::Int, tol::T) where {T}
        y = Vector{T}(undef, 2 * n + 1)
        w = Vector{T}(undef, 2 * n + 1)
        if n <= 0
            throw(DomainError(n, "the number of points must be > 0"))
        end
        h = convert(T, lambertW(Float64(2 * pi * n)) / n)
        for i = -n:n
            q = exp(-sinh(i * h)*pi)
            yi = 2 * q / (1 + q)
            y[n+i+1] = yi
            w[n+i+1] = yi * h * pi * cosh(i * h) / (1 + q)
            if isnan(w[n+i+1]) || w[n+i+1] < tol
                # overflow is expected to happen for large n
                # the correct way to address it is to set the weights to zero (ignore the points)
                w[n+i+1] = Base.zero(T)
                y[n+i+1] = one(T)
            end
        end
        return new{T}(h, y, w, tol)
    end
end

function integrate(q::TanhSinhQuadrature{T}, integrand)::T where {T}
    n = trunc(Int, (length(q.w) - 1) / 2)
    Il = Base.zero(T)
    Iprev = one(T)
    Ithreshold = eps(T)^2
    for i = n:-1:1
        yi = q.y[i]
        zi = (2 - yi) / yi
        if q.w[i] != 0
            fyi = integrand(zi)
            if fyi == 0 || isnan(fyi)
                break
            else
                uxi = 2 * fyi / yi^2
                Inew = q.w[i] * uxi
                # println(zi, " ",uxi)
                # println(i," ",zi, " ",Inew)
                Il += Inew
                if abs(Inew) + abs(Iprev) < Ithreshold
                    break
                else
                    Iprev = Inew
                end
            end
        end
    end
    Iu = Base.zero(T)
    Iprev = one(T)
    for i = n+1:2*n+1
        yi = q.y[i]
        zi = (2 - yi) / yi
        if q.w[i] != 0
            fyi = integrand(zi)
            if fyi == 0 || isnan(fyi)
                break
            else
                uxi = 2 * fyi / yi^2
                # println(zi, " ",uxi)
                Inew = q.w[i] * uxi
                Iu += Inew
                if abs(Inew) + abs(Iprev) < Ithreshold
                    break
                else
                    Iprev = Inew
                end
            end
        end
    end
    I = Iu + Il
end

function integrate(q::ModlobQuadrature{T}, integrand)::T where {T}
    #integrand with variable transform from -1 to 1
    @inline function integrand1(x::T)::T
        if x == 1
            return Base.zero(T)
        end
        u = (1 + x) / (1 - x)
        return integrand(u) * 2 / (1 - x)^2
    end
    I = modlob(integrand1, -one(T), one(T), tol = tolerance(q))
    return I
end


function integrate(q::DEQuadrature{T}, integrand)::T where {T}
    I, _ = quadde(integrand, 0, Inf, rtol = tolerance(q))
    return I
end



#Andersen and Lake "Robust High-Precision Option Pricing by Fourier Transforms: Contour Deformations and Double-Exponential Quadrature"
#Supports Heston & ShobelZhu
struct ALCharFuncPricer{T}
    cf::CharFunc
    quadrature::ALQuadrature{T}
end

function ALCharFuncPricer(cf:: Union{CharFunc{HestonParams{T}, CR},CharFunc{SchobelZhuParams{T},CR}}; n::Int=200) where {CR,T}
    return ALCharFuncPricer(cf, TanhSinhQuadrature(n,eps(T)))
end


function priceEuropean(
    p::ALCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    cf = p.cf
    m = model(cf)
    ω = log(forward / strike)
    r = m.ρ - m.σ * ω / (m.v0 + m.κ * m.θ * τ)
    piHigh = T(Base.pi)# pi(cf)
    angle = (r * ω < 0) ? piHigh / 12 * sign(ω) : Base.zero(T)
    α = optimalAlpha(p.cf, τ, ω)
    logmax = -log(floatmin(T)) - 8
    if α > 0
        α = min(α, logmax / abs(ω))
    elseif α < 0
        α = max(α, -logmax / abs(ω))
    end
    # println("alpha ", α, " om ", ω)
    tana = tan(angle)
    cc1 = oneim(cf)
    cc1a =  α * cc1
    cc1t = tana*cc1 +1
    nEval = 0
    @inline function integrand(x::T)::T
        hx = -cc1a + x * cc1t
        logChar = evaluateLogCharFunc(cf, hx - cc1, τ)
        # Qhx = exp(logChar) / (hx * (hx - cc1))
        # return real(exp(-x * tana * ω) * exp(cc1 * x * ω) * Qhx * (1 + cc1 * tana)) #exp may be merged into one
        v = exp(-tana * x * ω + logChar) / (hx * (hx - cc1)) * cc1t
        #real ( v* + exp(1im * x * ω))
        sx, cx = sincos(x * ω)
        nEval += 1
        return real(v) * cx - imag(v) * sx
    end
    I = Base.zero(T)
    if abs(ω * α) < logmax - sqrt(eps(T))
        I = integrate(p.quadrature, integrand)
        # println((α , I))
        I *= exp(α * ω) * forward / piHigh

    end
    #  println(" alpha ", α, " alphaom ", α * ω, " I ", I, " nEval ", nEval)
    R = Base.zero(T)
    if α <= 0
        R += forward
        if α == 0
            R -= forward / 2
        elseif α <= -1
            R -= strike
            if α == -1
                R += strike / 2
            end
        end
    end
    if isCall
        return (R - I) * discountDf
    else
        return ((R + strike - forward) - I) * discountDf
    end
end


@inline function kmp(
    cf::CharFunc{HestonParams{T},CR},
    τ::T,
    x::T,
)::Tuple{T,T} where {T,CR}
    p = model(cf)
    κ = p.κ
    σ = p.σ
    ρ = p.ρ
    squareTerm = (σ - 2 * κ * ρ)^2 + 4 * (κ^2 + x^2 / τ^2) * (1 - ρ^2)
    sqrtTerm = sqrt(max(zero(T), squareTerm))
    km = ((σ - 2 * κ * ρ) - sqrtTerm) / (2 * σ * (1 - ρ^2))
    kp = ((σ - 2 * κ * ρ) + sqrtTerm) / (2 * σ * (1 - ρ^2))
    return (km, kp)
end

@inline function evaluateβ(cf::CharFunc{HestonParams{T},CR}, k::T)::T where {T,CR}
    p = model(cf)
    β = p.κ - p.ρ * p.σ * k
    return β
end

@inline function evaluateD2(
    cf::CharFunc{HestonParams{T},CR},
    β::T,
    k::T,
)::T where {T,CR}
    p = model(cf)
    return D2 = β^2 - p.σ^2 * k * (k - 1)
end

function evaluateK(
    cf::Union{CharFunc{HestonParams{T},CR},CharFunc{SchobelZhuParams{T},CR}}, km0::T, kp0::T,
    τ::T,
    k::T,
)::T where {T,CR}
    β = evaluateβ(cf, k)
    D2 = evaluateD2(cf, β, k)
    if k >= km0 && k <= kp0
        D = sqrt(max(zero(T),D2))
        return D * τ - log((β - D) / (β + D))
    else
        D = sqrt(max(zero(T),-D2))
        if β > 0
            return D * τ - 2 * (Base.pi - atan(D / β))
        else
            return D * τ + 2 * (atan(D / β))
        end
    end
end


#roots of D(-ik) = -x^2 / τ^2
@inline function kmp(
    cf::CharFunc{SchobelZhuParams{T},CR},
    τ::T,
    x::T
)::Tuple{T,T} where {T,CR}
    p = model(cf)
    κ = p.κ
    σ = p.σ
    ρ = p.ρ
    squareTerm = (σ - 2 * κ * ρ)^2 + (4 * κ^2 + x^2 / τ^2) * (1 - ρ^2)
    sqrtTerm = sqrt(max(zero(T), squareTerm))
    km = ((σ - 2 * κ * ρ) - sqrtTerm) / (2 * σ * (1 - ρ^2))
    kp = ((σ - 2 * κ * ρ) + sqrtTerm) / (2 * σ * (1 - ρ^2))
    return (km, kp)
end

@inline function evaluateβ(cf::CharFunc{SchobelZhuParams{T},CR}, k::T)::T where {T,CR}
    p = model(cf)
    β = 2 * (p.κ - p.ρ * p.σ * k)
    return β
end

@inline function evaluateD2(
    cf::CharFunc{SchobelZhuParams{T},CR},
    β::T,
    k::T,
)::T where {T,CR}
    p = model(cf)
    return D2 = β^2 - 4 * p.σ^2 * k * (k - 1)
end

#ω is log(F/K)
function optimalAlpha(
    cf::Union{CharFunc{HestonParams{T},CR},CharFunc{SchobelZhuParams{T},CR}},
    τ::T,
    ω::T,
)::T where {T,CR}

    # θ = p.θ

    km0, kp0 = kmp(cf, τ, zero(T))
    #critial moments solution of M(k) == 1/D * log(beta-D / beta+D)  k in k-, k+ or ... (eqn 36)
    #solve in between brackets (proposition 2) via Brent.
    cc1 = oneim(cf)
    piHigh = const_pi(cf)


    km2pi, kp2pi = kmp(cf, τ, 2 * piHigh)
    ktol = T(1e-4)
    #we have kmin & kmax, \alpha = k-1
    @inline function alphaObj(α::T)::T
        # u = -(α + 1) * cc1
        logPhi = evaluateLogCharFuncAtImag(cf, -(α+1), τ)
        # println("alphaObj ",α," ",logPhi," ",evaluateLogCharFunc(cf, -(α + 1) * cc1, τ))
        if α >= -1 && α <= 0 #Andersen-Lake typical case, include denominator (payoff dependent)
            return real(logPhi) + α * ω - log(-α * (α + 1))
        else
            return real(logPhi) + α * ω - log(α * (α + 1))
        end
        #return real(logPhi) + α * ω
    end

    @inline function kObjective(x::T)::T
        value = evaluateK(cf, km0, kp0, τ, x)
        value
    end
    local αmin, αmax
    ismaxOk = true
    if ω < 0
        local kmax
        threshold = evaluateβ(cf, one(T))
        if threshold > 0
            if (kp2pi - 1) < 0.25
                kmax = one(T)
                ismaxOk = false ##too close to alpha=0
            else
                keps = ktol * min(1.0, kp2pi - kp0)
                kmax = find_zero(kObjective, (kp0 + keps, kp2pi - keps), A42())
            end
        elseif threshold < 0
            τCut = -2 / evaluateβ(cf, kp0)
            if τ < τCut
                kmpi, kppi = kmp(cf, τ, T(Base.pi))
                keps = ktol * min(one(T), kppi - kp0)
                kmax = find_zero(kObjective, (kp0 + keps, kppi - keps), A42())
            else
                keps = ktol * min(one(T), kp0 - 1)
                f0 = kObjective(one(T) + keps)
                f1 = kObjective(kp0 - keps)
                if sign(f1) == sign(f0) || kp0 - 1 < 0.25
                    kmax = one(T)
                    ismaxOk = false
                else
                    try
                        #kmax = find_zero(kObjective, (one(T)+keps, kp0 - keps), A42()) #sinhcosh formulation
                        kmax = find_zero(kObjective, (one(T) + keps, kp0 - keps), A42()) #log formulation
                    catch e
                        println(
                            km0,
                            " ",
                            kp0,
                            " ",
                            f0,
                            " ",
                            f1,
                            " ",
                            kObjective(1 + 10 * eps()),
                        )
                        rethrow(e)
                        # kmax = abs(f0) < abs(f1) ? 1 + keps : kp0
                        #    rethrow(e)
                        ismaxOk = false
                    end
                end
            end
        else #0
            kmpi, kppi = kmp(cf, τ, T(Base.pi))
            keps = ktol * min(1.0, kppi - kp0)
            kmax = find_zero(kObjective, (kp0 + keps, kppi - keps), A42())
        end
        αmin = Base.zero(T)
        αmax = kmax - 1
        if (αmax - αmin) < 0.25 #do not allow an alpha too close to zero
            ismaxOk = false
        end
    end
    if ω >= 0
        keps = ktol * (km0 - km2pi)
        #a low tolerance here turns out to be relatively important.
        kmin = find_zero(kObjective, (km2pi + keps, km0 - keps), A42())
        αmin = kmin - 1
        αmax = -one(T)
    end
    if !ismaxOk
        αmin = -one(T)
        αmax = Base.zero(T)
        #return -0.25
    end
    # println("αmin ",αmin, " αmax ",αmax)
    aeps = min(one(T), αmax - αmin) * ktol
    result = optimize(alphaObj, αmin + aeps, αmax - aeps, Brent(); rel_tol = T(1e-4))
    return minimizer(result)

end
