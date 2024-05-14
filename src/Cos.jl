export CosCharFuncPricer, makeCosCharFuncPricer, priceEuropean, priceDigital

struct CosCharFuncPricer{T}
    τ::T
    a::T
    b::T
    uk::Array{T,1}
    phi::Array{T,1}
    pi::T
end

#Gero Junike
using TaylorSeries
using ForwardDiff
using FastGaussQuadrature
function makeCosCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T;
    s=20, #number of derivatives to determine the number of terms
    relStrike=1.0,
    tol::T=1e-6
) where {T,CR,MAINT}
    t = Taylor1(Float64, 4)
    cft = evaluateCharFunc(cf, t, τ)
    phi0 = cft.coeffs[1]
    phi1 = cft.coeffs[2]
    mu = real(phi1 * -1im)
    phi2 = cft.coeffs[3] * 2
    phi3 = cft.coeffs[4] * 2 * 3
    phi4 = cft.coeffs[5] * 2 * 3 * 4
    phi4 = phi4 - 4im * mu * phi3 - 6 * mu^2 * phi2 + 4im * mu^3 * phi1 + mu^4 * phi0 #corresponds to (f(x)*exp(-i*x*mu))''''
    # phi4 = evaluateFourthDerivative(cf, τ, zero(T)); mu=0.0
    mu_n = abs(phi4) #4-th moment of log-returns. 
    l = (2 * relStrike * mu_n / tol)^(1 / 4) #Truncation range, Junike (2024, Eq. (3.10)).
    integrand = function (u)
        1 / (2 * pi) * abs(u)^(s + 1) * abs(evaluateCharFunc(cf, u, τ) * exp(-1im * u * mu))
    end
    if s <= 0
        c0 = cinf(model(cf), τ)
        lWinv = Float64(c0 / (tol^(3/2) * l))
        lW = lambertW(lWinv)
        m = ceil(Int, lW / Float64(c0) / Base.pi * Float64(2l))
        m = max(64, m)
    else
        # deResult = quadde(integrand, -Inf, Inf) #slow, gausshermite on 32 point does not work 
        al = ALTransformation()
        integrandT = function (z)
            integrand(transform(al, z)) * transformDerivative(al, z)
        end
        #println(integrandT(inverseTransform(al,1.0)))
        x, w = gausslegendre(128) #note: quadde on integrand does not work!
        f = u -> integrandT(u)
        boundDeriv = 2 * dot(w, f.(x))
        # z = @. (1:1024) * pi / 2l #if we want to reuse phi of cos method, unclear what N should be (circular problem), here 1024 is not enough
        # phiz = map(z -> evaluateCharFunc(cf, z, τ), z)
        # boundDeriv = 2*(z[2]-z[1])*sum(@.(abs(phiz .* exp(-1im*mu*z))*abs(z)^(s + 1)))/(2*pi)
        println("mu ",mu,"phi4 ",phi4, " mu_n ",mu_n, " l ",l," res ",boundDeriv)
        tmp = 2^(s + 5 / 2) * boundDeriv * l^(s + 2) * 12 * relStrike
        m = ceil(Int, (tmp / (s * pi^(s + 1) * tol))^(1 / s)) #Number of terms, Junike (2024, Sec. 6.1) 
    end
    return makeCosCharFuncPricer(cf, τ, m, -l, l)
end


#using cumulants rule
function makeCosCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T,
    m::Int,
    l::Int
) where {T,CR,MAINT}
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    # println("a ",a," b ",b)   
    return makeCosCharFuncPricer(cf, τ, m, a, b)
end

#On interval [a,b]
function makeCosCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T,
    m::Int,
    a::T,
    b::T
) where {T,CR,MAINT}
    p = model(cf)
    piHigh = const_pi(cf)
    z = @. (1:m) * piHigh / (b - a)
    phiz = map(z -> evaluateCharFunc(cf, z, τ), z)
    phi = @. real(phiz) * cos(-z * a) - imag(phiz) * sin(-z * a)
    uk = Vector{typeof(piHigh)}(undef, m)
    return CosCharFuncPricer(τ, a, b, uk, phi, piHigh)
end


function makeCosCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T,
    l::Int;
    tol::T=1e-8,
) where {T,CR,MAINT}
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    c0 = cinf(model(cf), τ)
    lWinv = Float64(2 * c0 / (tol * (b - a)))
    lW = lambertW(lWinv)
    m = ceil(Int, lW / Float64(c0) / Base.pi * Float64(b - a))
    m = max(32, m)
    # println("a ",a," b ",b, " m ",m, " lwi ",lWinv, " c0 ",c0)
    piHigh = const_pi(cf)
    z = @. (1:m) * piHigh / (b - a)
    phiz = map(z -> evaluateCharFunc(cf, z, τ), z)
    phi = @. real(phiz) * cos(-z * a) - imag(phiz) * sin(-z * a)
    uk = Vector{typeof(piHigh)}(undef, m)
    return CosCharFuncPricer(τ, a, b, uk, phi, piHigh)
end

#we adopt here the alternative formula of LeFloch "More Robust Pricing of European Options Based on Fourier Cosine Series Expansions"
function priceEuropean(
    p::CosCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    local pricePut
    x = log(forward / strike)
    if x >= -p.a && x >= p.b
        pricePut = 0
    elseif x <= p.a || x <= -p.b
        pricePut = discountDf * (strike - forward)
    else
        uk = p.uk
        a = p.a
        b = p.b
        ea = exp(a)
        f = forward
        piHigh = p.pi
        logStrike = log(strike / f)
        estrike = strike / f
        coeff = (-(estrike - ea) + estrike * (logStrike - a)) * 2 / (b - a)
        uk0 = coeff

        @inbounds for i = eachindex(uk)
            z = i * piHigh / (b - a)
            kPid = (logStrike - a) * z
            sk, ck = sincos(kPid)
            chi = (ck * estrike - 1 * ea + z * (sk * estrike)) / (1 + z^2)
            psi = sk / z
            coeff = (-chi + estrike * psi) * 2 / (b - a)
            uk[i] = coeff
        end
        sumPut = uk0 / 2
        @inbounds for k = eachindex(uk)
            sumPut += p.phi[k] * p.uk[k]
        end
        pricePut = discountDf * f * sumPut
    end
    if isCall
        return pricePut + discountDf * (forward - strike)
    end
    return pricePut
end


function priceDigital(
    p::CosCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    local pricePut
    x = log(forward / strike)
    if x >= -p.a && x >= p.b
        pricePut = 0
    elseif x <= p.a || x <= -p.b
        pricePut = discountDf
    else
        uk = p.uk
        a = p.a
        b = p.b
        ea = exp(a)
        f = forward
        piHigh = p.pi
        logStrike = log(strike / f)
        estrike = strike / f
        coeff = (logStrike - a) * 2 / (b - a)
        uk0 = coeff

        @inbounds for i = eachindex(uk)
            z = i * piHigh / (b - a)
            kPid = (logStrike - a) * z
            sk = sin(kPid)
            psi = sk / z
            coeff = psi * 2 / (b - a)
            uk[i] = coeff
        end
        sumPut = uk0 / 2
        @inbounds for k = eachindex(uk)
            sumPut += p.phi[k] * p.uk[k]
        end
        pricePut = discountDf * sumPut
    end
    if isCall
        return -pricePut + discountDf
    end
    return pricePut
end


using TaylorSeries
function computeCumulants(cf::CharFunc{MT,CR}, τ::T) where {MT,CR,T}
    t = Taylor1(T, 5)
    cft = evaluateLogCharFunc(cf, t, τ)
    return (imag(cft.coeffs[2]), -real(cft.coeffs[3]) * 2, real(cft.coeffs[5] * 2 * 3 * 4))
end

