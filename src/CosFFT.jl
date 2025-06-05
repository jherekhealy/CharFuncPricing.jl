struct CosFFTCharFuncPricer{MAINT,CR,T}
    cf::CharFunc{MAINT,CR}
    τ::T
    uk::Array{T,1}
    a::T
    b::T
end

function CosFFTCharFuncPricer(
    cf::CharFunc{MAINT,CR},
    τ::T,
    m::Int,
    l::Int
) where {T,CR,MAINT}
    c1, c2, c4 = computeCumulants(cf, τ)
    #println("c1 ",c1," c2 ",c2," c4 ",c4)
    c2 = abs(c2) + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))
    piHigh = const_pi(cf)
    # println("a ",a," b ",b)   
    return CosFFTCharFuncPricer(cf, τ, Vector{typeof(piHigh)}(undef, m), a, b)
end

using NFFT

function priceEuropean(
    p::CosFFTCharFuncPricer{MAINT,CR,T},
    isCalls::AbstractArray{Bool},
    strikes::AbstractArray{T},
    forward::T,
    τ::T,
    discountDf::T;
    isAlternative::Bool=true,
    useNFFT::Bool=true,
    reltol=1e-9
) where {MAINT,CR,T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    piHigh = const_pi(p.cf)
    ea = exp(p.a)
    m = length(p.uk)
    putPrices = if isAlternative
        i=1:m-1
        eta = @. piHigh * i / (p.b - p.a)
        uconst = @. 2/((p.b - p.a)*(1+eta^2))*ea
        upos = @. 2/((p.b - p.a)*(1+eta^2))*(-1 + 1im*eta)/2 + 2/((p.b - p.a)*eta)*(-1im/2)
        uneg = @. 2/((p.b - p.a)*(1+eta^2))*(-1 - 1im*eta)/2 + 2/((p.b - p.a)*eta)*(1im/2)
        phi = map(i -> evaluateCharFunc(p.cf, piHigh * i / (p.b - p.a), τ) * exp(-(p.a) * piHigh * i * 1im / (p.b - p.a)), 1:m-1)
        phi0 =  evaluateCharFunc(p.cf, zero(T), τ)
        constTerm = real(sum(phi .* uconst) + 2*ea/(p.b - p.a) *phi0 / 2)
        z = @. (log(strikes / forward)-p.a) / 2(p.b - p.a)   
        #FIXME make sure z is in [-0.5,0.5)
        fHat = if useNFFT
            f = vcat( zero(CR),reverse!(phi .* upos), zero(CR), phi .* uneg)
            real(nfft(z, f, reltol=reltol))     #, reltol=1e-16
        else
            map(x -> sum(real( phi[i] * (upos[i] * exp(i * piHigh * 1im * (2x)) + uneg[i]*exp(-i * piHigh * 1im * (2x))) ) for i = 1:m-1), z)
        end
        @. fHat = (forward*constTerm + strikes *real(phi0)* (z*2 - 1/(p.b-p.a)) +strikes * fHat )*discountDf
    else
        for i = 1:m
            etak = piHigh * i / (p.b - p.a)
            sk, ck = sincos(etak * p.a)
            p.uk[i] = ((ea + etak * sk - ck) / (etak^2 + 1) - sk / etak) * 2 / (p.b - p.a)
        end
        uk0 = (ea - 1 - p.a) * 2 / (p.b - p.a)
        #sumPut = real(evaluateCharFunc(p.cf, zero(T), τ)) * uk0 / 2
        # x in [a,b] means pi*(x-a)/(b-a) in [0,pi]
        # or 2*pi*n*z with z in [0,1/2] and n in 1,N-1
        # if we set -b=a, we have (x-a)/b-a = x/2b x-a/b-a =  x/2b + 1/2 = z+1/2 with z in [-1/2,1/2] n in [1,N-1] or -2*k*(z+1/2) with k in [1/2, N-1/2]
        f0 = evaluateCharFunc(p.cf, zero(T), τ) * uk0 / 2
        phi = map(i -> evaluateCharFunc(p.cf, piHigh * i / (p.b - p.a), τ) * exp(-(p.a) * piHigh * i * 1im / (p.b - p.a)), 1:m-1)

        z = @. log(strikes / forward) / 2(p.b - p.a)
        #FIXME make sure z is in [-0.5,0.5)
        fHat = if useNFFT
            f = vcat(zeros(CR, m), f0, phi .* p.uk[1:m-1])
            real(nfft(z, f, reltol=reltol))     #, reltol=1e-16
        else
            map(x -> real(f0) + sum(real(p.uk[i] * phi[i] * exp(i * piHigh * 1im * (-2x))) for i = 1:m-1), z)
        end
        @. fHat *= discountDf * strikes
        fHat
    end

    return map((isCall, strike, pricePut) -> ifelse(isCall, pricePut + discountDf * (forward - strike), pricePut), isCalls, strikes, putPrices)
end

