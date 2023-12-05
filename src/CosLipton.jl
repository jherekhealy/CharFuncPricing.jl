export CosLiptonCharFuncPricer, makeCosLiptonCharFuncPricer, priceEuropean

struct CosLiptonCharFuncPricer{T}
    τ::T
    a::T
    b::T
    uk::AbstractArray{T}
    phi::AbstractArray{Complex{T}}
    pi::T
end


function makeCosLiptonCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T,
    m::Int,
    l::Int,
) where {T,CR,MAINT}
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    # println("a ",a," b ",b)
    piHigh = const_pi(cf)
    z = @. (0:m) * piHigh / (b - a)
    phiz = map(z -> evaluateCharFunc(cf, z-0.5im, τ), z)
    uk = Vector{typeof(piHigh)}(undef, m)
    return CosLiptonCharFuncPricer(τ, a, b, uk, phiz, piHigh)
end


#we adopt here the alternative formula of LeFloch "More Robust Pricing of European Options Based on Fourier Cosine Series Expansions"
function priceEuropean(
    p::CosLiptonCharFuncPricer{T},
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
        cha = cosh(a/2)
        if isCall
            cha = cosh(b/2)
        end
        piHigh = p.pi        
        uk0 = (cha-1)*4
        # uk0 = exp(a)-1-a
        # uk0 =  (exp(b) -1)-b
        sumPut = uk0/2*real(p.phi[1])*2 / (b - a)

        @inbounds for i = 1:length(uk)
            z = i * piHigh / (b - a)
            kPid = ( a) * z
            ck = cos(kPid)
            sign = 1
            if isCall && i % 2 == 1
                sign = -1
            end
            uk[i] = (cha*sign - ck) / (0.25 + z^2)
            # sk,ck = sincos(kPid)
            # uk[i] = (exp(a)+z*sk-ck)/(1+z*z)-sk/z
            
            # uk[i] = (exp(b)*(-1)^i -ck +z*sk)/(1+z*z)-sk/z
            phi = real(p.phi[i+1]) * cos(-z * (a-x)) - imag(p.phi[i+1]) * sin(-z * (a-x))
            sumPut += phi * uk[i]*2 / (b - a)
        end
        # sumPut *= 
        pricePut = discountDf  * sqrt(forward* strike) * sumPut
        # pricePut = discountDf  *  strike * sumPut
    end
    # if isCall
    #     return pricePut + discountDf * (forward - strike)
    # end
    return pricePut
end
