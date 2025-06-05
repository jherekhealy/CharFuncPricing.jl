using FCCQuad
#TODO add alpha, and price strikes together, no need for cache then.

struct FCCACharFuncPricer{MAINT,T,CR}
    cf::CharFunc{MAINT,CR}
    τ::T
    b::T
    vControl::T
    qTol::T
    α::T
    function FCCACharFuncPricer(
        p::CharFunc{MAINT,CR},
        τ::T;
        b::T=Base.zero(T),
        α = -0.5,
        tTol::T=sqrt(eps(T)),
        qTol::T=sqrt(eps(T)),
    ) where {MAINT,CR,T}
        if b == 0
            b = computeTruncation(p, τ, tTol)
        end
        return new{MAINT, T,CR}(p,τ, b, getControlVariance(p), qTol, α)
    end
end

Base.broadcastable(p::FCCACharFuncPricer) = Ref(p)

function priceEuropean(
    p::FCCACharFuncPricer{MAINT, T,CR},
    isCalls::AbstractArray{Bool},
    strikes::AbstractArray{T},
    forward::T,
    τ::T,
    discountDf::T,
) where {MAINT, CR, T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    freqs = @. -log(strikes / forward) #-k
    iPure = oneim(p.cf)
    @inline function integrand(x::T)::CR where {T}
        u = x - (p.α+1)*iPure
        phi = evaluateCharFunc(p.cf, u, τ)
        denom = (x * iPure +p.α)*(x*iPure + p.α + 1)
        phi /= denom
        #println("integrand x=", x, " phi=", phi, " denom=", denom)
        return phi
    end
    
    resultAndError, neval = fccquad(integrand,  x -> one(T), freqs, xmin=zero(T), xmax=p.b, method=:degree, abstol=p.qTol)
    #println("neval=", neval, " resultAndError=", resultAndError)
    callPrices = map((freq,integral) -> forward* exp(freq*p.α) / pi * real(integral), freqs, resultAndError[1,:])
    if p.vControl != 0
        #use control variate
        variance = p.vControl * p.τ
        sqrtVar = sqrt(variance)
        d1 = log(forward / strike) / sqrtVar + sqrtVar / 2
        d2 = d1 - sqrtVar
        nd1 = normcdf(d1)
        nd2 = normcdf(d2)
        price = discountDf * (forward * nd1 - strike * nd2)
        callPrice += price
    else
        #FIXME case where alpha >0
        if p.α <= 0
            @. callPrices += forward
        end
        if p.α <= -1
            @. callPrices -= strikes
        else
            if p.α == 0
                @. callPrices -= forward/2
            elseif p.α == -1
                @. callPrices += strikes/2
            end
        end
    end
    return map( (isCall,callPrice,strike) -> ifelse(isCall, callPrice*discountDf, (callPrice - (forward - strike)) * discountDf), isCalls, callPrices, strikes)
end

struct FCCCharFuncPricer{T,CR}
    τ::T
    b::T
    vControl::T
    qTol::T
    cache::Vector{CR}
    function FCCCharFuncPricer(
        p::CharFunc{MAINT,CR},
        τ::T;
        b::T=Base.zero(T),
        tTol::T=sqrt(eps(T)),
        qTol::T=sqrt(eps(T)),
    ) where {MAINT,CR,T}
        if b == 0
            b = computeTruncation(p, τ, tTol)
        end
        iPure = oneim(p)
        cache = Vector{CR}()
        @inline function integrand(x::T)::CR where {T}
            u = x - 0.5 * iPure
            phi = evaluateCharFunc(p, u, τ)
            denom = 0.25 + x^2
            phi /= denom
            #println("integrand x=", x, " phi=", phi, " denom=", denom)
            push!(cache, phi)
            return phi
        end
        fccquad(integrand, x -> one(T), [zero(T)], xmin=zero(T), xmax=b, method=:degree, abstol=qTol)
        return new{T,CR}(τ, b, getControlVariance(p), qTol, cache)
    end
end

Base.broadcastable(p::FCCCharFuncPricer) = Ref(p)

function priceEuropean(
    p::FCCCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    freq = -log(strike / forward) #-k
    index = 0
    resultAndError, neval = fccquad(function (x)
            index += 1
            return p.cache[index]
        end, x -> one(T), [freq], xmin=zero(T), xmax=p.b, method=:degree, abstol=p.qTol)
    integral = real(resultAndError[1])
    callPrice = -discountDf * sqrt(strike * forward) / pi * integral
    if p.vControl != 0
        #use control variate
        variance = p.vControl * p.τ
        sqrtVar = sqrt(variance)
        d1 = log(forward / strike) / sqrtVar + sqrtVar / 2
        d2 = d1 - sqrtVar
        nd1 = normcdf(d1)
        nd2 = normcdf(d2)
        price = discountDf * (forward * nd1 - strike * nd2)
        callPrice += price
    else
        callPrice += discountDf * forward
    end
    if isCall
        return callPrice
    end
    putPrice = callPrice - (forward - strike) * discountDf
    return putPrice
end