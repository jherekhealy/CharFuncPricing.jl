#Cos pricing based on integrating on the -i line. Calls can be priced, and Puts must use put-call parity.
struct CosGPCharFuncPricer{T}
    τ::T
    a::T
    b::T
    uk::Array{T,1}
    phi::Array{T,1}
    phiF::Array{T,1}
    pi::T
end


#using cumulants rule
function makeCosGPCharFuncPricer(
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
    return makeCosGPCharFuncPricer(cf, τ, m, a, b)
end

#On interval [a,b]
function makeCosGPCharFuncPricer(
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
    phii = evaluateCharFunc(cf, zero(T)-1im, τ)
    phizF = map(z -> evaluateCharFunc(cf, z-1im, τ)/phii, z)
    phi = @. real(phiz) * cos(-z * a) - imag(phiz) * sin(-z * a)
    phiF = @. real(phizF) * cos(-z * a) - imag(phizF) * sin(-z * a)
    
    uk = Vector{typeof(piHigh)}(undef, m)
    return CosGPCharFuncPricer(τ, a, b, uk, phi, phiF, piHigh)
end



function priceEuropean(
    p::CosGPCharFuncPricer{T},
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
         sumPutF = uk0 / 2
        @inbounds for k = eachindex(uk)
            sumPutF += p.phiF[k] * p.uk[k]
        end
        pricePut = discountDf * (strike*sumPut - forward*sumPutF)
    end
    if isCall
        return (pricePut +forward-strike)* discountDf
    end
    return pricePut
end



