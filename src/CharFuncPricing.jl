module CharFuncPricing
include("Heston.jl")

export CosCharFuncPricer, makeCosCharFuncPricer, priceEuropean

struct CosCharFuncPricer{T,CT}
    τ::T
    a::T
    b::T
    uk::Vector{T}
    phi::Vector{CT}
    pi::T
end


function makeCosCharFuncPricer(CC, R, pi::T, p::HestonParams{T}, τ::T, m::Int, l::Int) where {T}
    c1, c2, c4 = computeCumulants(p, τ)
    c2 = c2 + sqrt(abs(c4))
    a = c1 - l * sqrt(abs(c2))
    b = c1 + l * sqrt(abs(c2))

    phi = zeros(CC, m)
    #phi[1] = evaluateCharFunc(p,0.0,τ)
    #V = collect(range(vmin,stop=vmax,length=L))
    for i = 1:m
        z = i * pi / (b - a)
        phi[i] = evaluateCharFunc(CC, p, z, τ)
    end
    uk = zeros(R, m)
    return CosCharFuncPricer(τ, a, b, uk, phi, pi)
end

function priceEuropean(p::CosCharFuncPricer{T,CT}, isCall::Bool, strike::T, forward::T, discountDf::T) where {T,CT}
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
        x = 0 #log(forward / f)

        phi0 = 1 #real(cfi.phi[0])
        sumPut = (phi0) * uk0 / 2
        for k = 1:length(uk)
            z = k * pi / (p.b - p.a)
            sz = sin(z * (x - p.a))
            cz = cos(z * (x - p.a))
            rephi = real(p.phi[k]) * cz - imag(p.phi[k]) * sz
            sumPut += rephi * p.uk[k]
        end
        pricePut = discountDf * f * sumPut
    end
    if isCall
        return pricePut + discountDf * (forward - strike)
    end
    return pricePut
end

end # module
