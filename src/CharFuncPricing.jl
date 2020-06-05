module CharFuncPricing
include("Heston.jl")

export CosCharFuncPricer, makeCosCharFuncPricer, priceEuropean

struct CosCharFuncPricer{T}
    τ::T
    a::T
    b::T
    uk::Vector{T}
    phi::Vector{T}
    pi::T
end


function makeCosCharFuncPricer(CC, R, pi::T, p::HestonParams{T}, τ::T, m::Int, l::Int) where {T}
    c1, c2, c4 = computeCumulants(p, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    phi = zeros(R, m)
    for i = 1:m
        z = i * pi / (b - a)
        phiz = evaluateCharFunc(CC, p, z, τ)
        sz = sin(-z * a)
        cz = cos(-z * a)
        phi[i] = real(phiz) * cz - imag(phiz) * sz
    end
    uk = zeros(R, m)
    return CosCharFuncPricer(τ, a, b, uk, phi, pi)
end

#we adopt here the alternative formula of LeFloch "More Robust Pricing of European Options Based on Fourier Cosine Series Expansions"
function priceEuropean(p::CosCharFuncPricer{T}, isCall::Bool, strike::T, forward::T, discountDf::T) where {T,CT}
    pricePut = 0
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
        pi = p.pi
        logStrike = log(strike / f)
        estrike = strike / f
        coeff = 2 / (b - a) * (-(estrike - ea) + estrike * (logStrike - a))
        uk0 = coeff

        for i = 1:length(uk)
            z = i * pi / (b - a)
            kPid = (logStrike - a) * z
            sk = sin(kPid)
            ck = cos(kPid)
            chi = 1 / (1 + z * z) * (ck * estrike - 1 * ea + z * (sk * estrike))
            psi = sk / z
            coeff = 2 / (b - a) * (-chi + estrike * psi)
            uk[i] = coeff
        end
        phi0 = 1 #real(cfi.phi[0])
        sumPut = (phi0) * uk0 / 2
        for k = 1:length(uk)
            sumPut += p.phi[k] * p.uk[k]
        end
        pricePut = discountDf * f * sumPut
    end
    if isCall
        return pricePut + discountDf * (forward - strike)
    end
    return pricePut
end

end # module
