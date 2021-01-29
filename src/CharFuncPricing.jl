module CharFuncPricing
include("lambertw.jl")
include("Heston.jl")
include("SchobelZhu.jl")
export CosCharFuncPricer, makeCosCharFuncPricer, priceEuropean

struct CosCharFuncPricer{T}
    τ::T
    a::T
    b::T
    uk::Array{T,1}
    phi::Array{T,1}
    pi::T
    field
end


function makeCosCharFuncPricer(
    CC,
    R,
    pi::T,
    p::Union{HestonParams{T},SchobelZhuParams{T}},
    τ::T,
    m::Int,
    l::Int,
) where {T}
    c1, c2, c4 = computeCumulants(p, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    # println("a ",a," b ",b)
    z =  @. (1:m) * pi / (b - a)
    phiz = map(z -> evaluateCharFunc(CC, p, z, τ), z)
    phi = @. real(phiz) * cos(-z * a) - imag(phiz) * sin(-z * a)
    uk = zeros(R, m)
    return CosCharFuncPricer(τ, a, b, uk, phi, pi,R)
end

#we adopt here the alternative formula of LeFloch "More Robust Pricing of European Options Based on Fourier Cosine Series Expansions"
function priceEuropean(
    p::CosCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    discountDf::T,
) where {T}
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
        pi = p.pi
        logStrike = log(strike / f)
        estrike = strike / f
        coeff = 2 / (b - a) * (-(estrike - ea) + estrike * (logStrike - a))
        uk0 = coeff

        @inbounds for i = 1:length(uk)
            z = i * pi / (b - a)
            kPid = (logStrike - a) * z
            sk, ck = sincos(kPid)
            chi = (ck * estrike - 1 * ea + z * (sk * estrike)) / (1 + z^2)
            psi = sk / z
            coeff = 2 / (b - a) * (-chi + estrike * psi)
            uk[i] = coeff
        end
        phi0 = 1
        sumPut = phi0 * uk0 / 2
        @inbounds @simd for k = 1:length(uk)
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
