export CosAlphaCharFuncPricer, makeCosAlphaCharFuncPricer, priceEuropean

struct CosAlphaCharFuncPricer{T}
    τ::T
    a::T
    b::T
    α::T
    uk::AbstractArray{T}
    phi::AbstractArray{Complex{T}}
    pi::T
end


function makeCosAlphaCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T,
    m::Int,
    l::Int,
    α::T
) where {T,CR,MAINT}
    p = model(cf)
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    # println("a ",a," b ",b)
    piHigh = const_pi(cf)
    z = @. (0:m) * piHigh / (b - a)
    phiz = map(z -> evaluateCharFunc(cf, z-Complex(0,α+1), τ), z)
    uk = Vector{typeof(piHigh)}(undef, m)
    return CosAlphaCharFuncPricer(τ, a, b, α, uk, phiz, piHigh)
end


#we adopt here the alternative formula of LeFloch "More Robust Pricing of European Options Based on Fourier Cosine Series Expansions"
function priceEuropean(
    p::CosAlphaCharFuncPricer{T},
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
        uk = p.uk
        α = p.α
        a = p.a 
        b = p.b 
        piHigh = p.pi        
        sumPut = zero(T)
    if isCall 
        eαp1 = exp(-α*b)
        eα = exp(-(α+1)*b)
        # uk0 = exp(a)-1-a
        # uk0 =  (exp(b) -1)-b
        χ = (-(α+1)*eα + (α+1)) / ((α+1)^2)
        χp1 = if α == zero(T)
            b   #(1-alpha*b) * -(alpha) + alpha  / alpha2 
        else 
            (-(α)*eαp1 + (α)) / ((α)^2)  
        end 
        uk0 = χp1-χ
        sumPut = uk0*real(p.phi[1]) / (b - a)
        @inbounds for i = 1:length(uk)
            z = i * piHigh / (b - a)
            kPid = ( a) * z
            sk, ck = sincos(kPid)
            χ = (-(α+1)*(-1)^i*eα + (α+1)*ck + z*sk) / ((α+1)^2 + z^2)
            χp1 = (-(α)*(-1)^i*eαp1 + (α)*ck + z*sk) / ((α)^2 + z^2)
            uk[i] = χp1-χ
            # sk,ck = sincos(kPid)
            # uk[i] = (exp(a)+z*sk-ck)/(1+z*z)-sk/z
            
            # uk[i] = (exp(b)*(-1)^i -ck +z*sk)/(1+z*z)-sk/z
            phi = real(p.phi[i+1]) * cos(-z * (a-x)) - imag(p.phi[i+1]) * sin(-z * (a-x))
            sumPut += phi * uk[i]*2 / (b - a)
        end
    else 
        eαp1 = exp(-α*a)
        eα = exp(-(α+1)*a)
        #println("eα ",eα, " ep1 ", eαp1, " eax ",exp((α+1)*x) )
        # uk0 = exp(a)-1-a
        # uk0 =  (exp(b) -1)-b
        χ = (-(α+1)*eα + (α+1)) / ((α+1)^2)
        χp1 = (-(α)*eαp1 + (α)) / ((α)^2)
        uk0 = χp1-χ
        sumPut = uk0*real(p.phi[1]) / (b - a)
        @inbounds for i = 1:length(uk)
            z = i * piHigh / (b - a)
            kPid = ( a) * z
            sk, ck = sincos(kPid)
            χ = (-(α+1)*eα + (α+1)*ck + z*sk) / ((α+1)^2 + z^2)
            χp1 = (-(α)*eαp1 + (α)*ck + z*sk) / ((α)^2 + z^2)
            uk[i] = χp1-χ
            # sk,ck = sincos(kPid)
            # uk[i] = (exp(a)+z*sk-ck)/(1+z*z)-sk/z
            
            # uk[i] = (exp(b)*(-1)^i -ck +z*sk)/(1+z*z)-sk/z
            phi = real(p.phi[i+1]) * cos(-z * (a-x)) - imag(p.phi[i+1]) * sin(-z * (a-x))
            sumPut += phi * uk[i]*2 / (b - a)
        end
    end
        # sumPut *= 
        price = discountDf  * strike* exp((α+1)*x) * sumPut
        # pricePut = discountDf  *  strike * sumPut
    
end
