export CosCharFuncPricer, makeCosCharFuncPricer, priceEuropean, priceDigital

struct CosCharFuncPricer{T}
    τ::T
    a::T
    b::T
    uk::Array{T,1}
    phi::Array{T,1}
    pi::T
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
        pricePut = discountDf * f * sumPut
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

