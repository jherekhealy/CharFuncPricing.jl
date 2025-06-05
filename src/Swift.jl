using FFTW, LinearAlgebra
export SwiftCharFuncPricer

struct SwiftCharFuncPricer{T}
    τ::T
    c::T
    m::Int
    j::Int
    barj::Int
    k1::Int
    k2::Int
    uk::Array{T,1}
    cl::Array{T,1}
end

asymptoticLogOrder(model::Union{HestonParams{T},SchobelZhuParams{T},SVCJParams{T}}) where {T} = one(T)


function findSwiftScaling(cf::CharFunc{MAINT,CR}, τ::T; mGuess=2, tol=sqrt(eps(T)), isInverseTime=false, mMax=20) where {T,CR,MAINT}
    #find m, given nu, which is model dependent. Maree uses inverseTime factor, but Leitao, Ortiz Garcia 2018 does not 
    ν = asymptoticLogOrder(model(cf))
    epsilon = tol
    mGuess -= 1
    while epsilon >= tol && mGuess < mMax
        mGuess += 1
        c = 2^mGuess * pi
        epsilon = (c)^(1 - ν) / (2 * pi * ν) * (abs(evaluateCharFunc(cf, -c, τ)) + abs(evaluateCharFunc(cf, c, τ)))
        if isInverseTime
            epsilon /= τ
        end
        # println(mGuess, " ",epsilon)
    end
    #in theory this may depend on x = ln(S(T)/K)
    return mGuess, epsilon
end

#TODO constructor with m optimized from init guess. Useful for more challenging params. Also l gives a first rough guess of accuracy.
SwiftCharFuncPricer(cf::CharFunc{MAINT,CR}, τ::T, m::Int, l::Number; tol=sqrt(eps(T))) where {T,CR,MAINT} = makeSwiftCharFuncPricer(cf, τ,m,l, tol=tol)
    
function makeSwiftCharFuncPricerMaree(cf::CharFunc{MAINT,CR}, τ::T, m::Int, l::Number; tol=sqrt(eps(T)), densityCalculator=TrapezoidalDensityCalculator()) where {T,CR,MAINT}
    #keep c, m, increase kappa in line with Maree
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = abs(c2) + sqrt(abs(c4))
    c = abs(c1) + l * sqrt(abs(c2)) #should also be larger that log(strikemax/f)
    κ = ceil(Int, 2^m * c)
    tol1 = tol
    twom2 = 2^(m / 2)
    J = 0
    local cl
    local phiz = zeros(CR, 0)
    while tol1 >= tol && J < 2^16
        if J != 0
            κ = ceil(Int, κ * 2) #what if κ is too large? must start with smallish estimate; we may cap it from J=16 just in case
        end
        ln2J = ceil(Int, log2(pi * κ))
        J = 2^ln2J
        if length(phiz) == J + 1
            #do nothing, reuse phiz
        elseif length(phiz) > 0
            phiold = phiz
            phiz = zeros(CR, J + 1)
            phiz[1:2:end] .= phiold
        else
            phiz = zeros(CR, J + 1)
        end
        cl = calculate(densityCalculator, cf, τ, m, κ, J, phiz)
        tol1 = abs(1 - sum(cl) / twom2)
        println(ln2J, " tol1=", tol1) #in master thesis, m is incremented and kappa deduced from ceil. 
    end
    return SwiftCharFuncPricer(τ, c, m, J, J, κ, κ, zeros(T, length(cl)), cl)
end

function makeSwiftCharFuncPricer(cf::CharFunc{MAINT,CR}, τ::T, m::Int, l::Number; tol=sqrt(eps(T)), densityCalculator=TrapezoidalDensityCalculator(), Jfactor=1) where {T,CR,MAINT}
    #increase c, keep m constant, deduce j from c. More in line with Leitao ORtiz.
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = abs(c2) + sqrt(abs(c4))
    c = abs(c1) + l * sqrt(abs(c2)) #should also be larger that log(strikemax/f)
    κ = ceil(Int, 2^m * c)
    twom2 = 2^(m / 2)
    tol1 = tol
    J = 0
    local cl
    local κ
    local phiz = zeros(CR, 0)
    while tol1 >= tol && J < 2^22
        if J != 0
            c *= sqrt(sqrt(2))
        end
        κ = ceil(Int, 2^m * c)
        ln2J = ceil(Int, log2(κ))    # the paper advise log2(pi*κ) because of particular series expansion. log2(κ) is sufficient for FFT.    
        J = 2^ln2J
        lenPhi = Jfactor * J + 1
        if length(phiz) == lenPhi
            #do nothing, reuse phiz
        elseif length(phiz) > 0
            phiold = phiz
            phiz = zeros(CR, lenPhi)
            phiz[1:2:end] .= phiold
        else
            phiz = zeros(CR, lenPhi)
        end
        cl = calculate(densityCalculator, cf, τ, m, κ, lenPhi - 1, phiz)
        tol1 = abs(1 - (sum(cl[1:end]) - cl[κ] / 2 - cl[κ+1] / 2) / twom2)
        # tol1 = abs(1 - (sum(cl[1:end])) / twom2)
        # println(m, " ", κ, " ", ln2J, " ", tol1)
    end    
    # c = (κ-1)/2^m #this is what Ortiz-Oosterlee would use for a put.
    return SwiftCharFuncPricer(τ, c, m, J, J, κ, κ, zeros(T, length(cl)), cl)
  
    #|y|=ln(S(T)/K) must be larger than c when interval is -c, c, 
  
end


function makeSwiftCharFuncPricerRomo(cf::CharFunc{MAINT,CR}, τ::T, m::Int, l::Number; tol=sqrt(eps(T)), densityCalculator=TrapezoidalDensityCalculator()) where {T,CR,MAINT}
    #increase c, keep m constant, deduce j from c. More in line with Leitao ORtiz.
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = abs(c2) + sqrt(abs(c4))
    c = abs(c1) + l * sqrt(abs(c2)) #should also be larger that log(strikemax/f)
    tol1 = tol
    J = 0
    local cl
    local κ
    local phiz = zeros(CR, 0)
    while tol1 >= tol && J < 2^20
        if J != 0
            m += 1
        end
        κ = ceil(Int, 2^m * c)
        twom2 = 2^(m / 2)
        ln2J = ceil(Int, log2(pi * κ))
        J = 2^ln2J

        phiz = zeros(CR, J + 1)
        cl = calculate(densityCalculator, cf, τ, m, κ, J, phiz)
        tol1 = abs(1 - (sum(cl[1:end]) - cl[κ] / 2 - cl[κ+1] / 2) / twom2)
        # tol1 = abs(1 - (sum(cl[1:end])) / twom2)
        println(m, " ", κ, " ", ln2J, " ", tol1) #in master thesis, m is incremented and kappa deduced from ceil. 
    end
    return SwiftCharFuncPricer(τ, c, m, J, J, κ, κ, zeros(T, length(cl)), cl)
    #

    #|y|=ln(S(T)/K) must be larger than c when interval is -c, c, 
    #this means kappa > 2^m c and we also must have J >= pi*kappa    logJ = log(pi*kappa)
end
abstract type DensityCoeffCalculator end

struct VietaDensityCalculator <: DensityCoeffCalculator
end

function calculate(calc::VietaDensityCalculator, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, J::Int, phiz::AbstractArray{CR}) where {T,MAINT,CR}
    piHigh = const_pi(cf)
    ω = @. (((1:J) - 1 / 2) * piHigh / J)
    z = ω * 2^m
    N = 2J

    for (j, zj) = enumerate(z)
        phiz[j] = evaluateCharFunc(cf, -zj, τ)
    end
    clz = zeros(CR, N)
    clz[1:length(phiz)] .= phiz
    fft!(clz)
    factor = 2^(m / 2) / J
    cl = zeros(T, 2κ)
    for k = 1:κ
        cl[2κ-k+1] = real(clz[k] * exp(-1im * piHigh * (k - 1) / N)) * factor
    end
    for k = -κ+1:0
        cl[-k+1] = real(clz[N+k] * exp(-1im * piHigh * (k - 1) / N)) * factor
    end

    return cl
end

struct SimpsonDensityCalculator <: DensityCoeffCalculator
end

function calculate(calc::SimpsonDensityCalculator, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, Jd::Int, phiz::AbstractArray{CR}) where {T,MAINT,CR}
    # J = Jd >> 1 #if we halve, the error may be unexpectedly large.
    clT = calculate(TrapezoidalDensityCalculator(), cf, τ, m, κ, J, zeros(CR, J+ 1))
    clV = calculate(VietaDensityCalculator(), cf, τ, m, κ, J, zeros(CR, J))
    cl = zeros(T, 2κ)
    for k = 1:2κ
        cl[k] = (clT[k] + 2 * clV[k]) / 3
    end
    return cl
end

struct SimpsonSlowDensityCalculator <: DensityCoeffCalculator
end

function calculate(calc::SimpsonSlowDensityCalculator, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, J::Int, phiz::AbstractArray{CR}) where {T,MAINT,CR}
    piHigh = const_pi(cf)
    ω = @. (((1:J)) * piHigh / J)
    z = ω * 2^m
    phiz0 = real(exp(evaluateLogCharFunc(cf, zero(T), τ)))
    for (j, zj) = enumerate(z)
        phiz[j] = evaluateCharFunc(cf, -zj, τ) #equivalently w could evaluate at zj and multiply by exp(-ikw) instead of exp(ikw)
    end
    cl = zeros(T, 2 * κ)
    factor = 2^(m / 2) / 3J
    for k = 1:κ
        cl[k] = (phiz0 + real(phiz[J] *exp( 1im * k * piHigh)) + 4sum([real(phizj * exp( 1im * k * ωj)) for (phizj, ωj) = zip(phiz[1:2:J-1], ω[1:2:J-1])]) + 2sum([real(phizj*exp( 1im * k * ωj)) for (phizj, ωj) = zip(phiz[2:2:J-2], ω[2:2:J-2])])) * factor
    end
    for k = -κ+1:0
        cl[k+2κ] = (phiz0 + real(phiz[J]*exp( 1im * k * piHigh)) + 4sum([real(phizj * exp( 1im * k * ωj)) for (phizj, ωj) = zip(phiz[1:2:J-1], ω[1:2:J-1])]) + 2sum([real(phizj*exp( 1im * k * ωj)) for (phizj, ωj) = zip(phiz[2:2:J-2], ω[2:2:J-2])])) * factor
    end

    return cl
end
struct TrapezoidalDensityCalculator <: DensityCoeffCalculator
end

function calculate(calc::TrapezoidalDensityCalculator, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, J::Int, phiz::AbstractArray{CR}) where {T,MAINT,CR}
    piHigh = const_pi(cf)
    ω = @. (((0:J)) * piHigh / J)
    z = ω * 2^m
    N = 2J
    for (j, zj) = enumerate(z)
        if phiz[j] == zero(CR)
            phiz[j] = evaluateCharFunc(cf, -zj, τ)
        end
    end
    clz = zeros(CR, N)
    clz[1:length(phiz)] .= phiz
    clz[1] /= 2
    clz[J+1] /= 2 #necessary as trapezoidal is on parseval identity, not an expansion.
    clz = ifft!(clz)
    factor = 2^(m / 2) * 2 #ifft has /N factor and N/J = 2
    cl = zeros(T, 2κ)
    for k = 1:κ
        cl[k] = real(clz[k+1]) * factor
    end
    for k = -κ+1:-1
        cl[2κ+k] = real(clz[N+k+1]) * factor
    end
    cl[2κ] = real(clz[1]) * factor
    return cl
end
struct TrapezoidalSlowDensityCalculator <: DensityCoeffCalculator
end

function calculate(calc::TrapezoidalSlowDensityCalculator, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, J::Int, phiz::AbstractArray{CR}) where {T,MAINT,CR}
    piHigh = const_pi(cf)
    ω = @. (((1:J)) * piHigh / J)
    z = ω * 2^m
    N = 2J
    phiz0 = real(evaluateCharFunc(cf, zero(T), τ))
    for (j, zj) = enumerate(z)
        phiz[j] = evaluateCharFunc(cf, -zj, τ) #equivalently w could evaluate at zj and multiply by exp(-ikw) instead of exp(ikw)
    end
    cl = zeros(T, 2 * κ)
    factor = 2^(m / 2) / J
    for k = 1:κ
        cl[k] = (phiz0 / 2 - real(phiz[J]*exp( 1im * k * piHigh)) / 2 + sum([real(phizj *exp(1im * k * ωj)) for (phizj, ωj) = zip(phiz, ω)])) * factor
    end
    for k = -κ+1:0
        cl[k+2κ] = (phiz0 / 2 - real(phiz[J]*exp( 1im * k * piHigh)) / 2 + sum([real(phizj*exp( 1im * k * ωj)) for (phizj, ωj) = zip(phiz, ω)])) * factor
    end

    return cl
end

struct EulerMaclaurinDensityCalculator{Order} <: DensityCoeffCalculator
end

function calculate(calc::EulerMaclaurinDensityCalculator{Order}, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, J::Int, phiz::AbstractArray{CR}) where {T,MAINT,CR,Order}
    cl = calculate(TrapezoidalDensityCalculator(), cf, τ, m, κ, J, phiz)
    if Order == 1
        piHigh = const_pi(cf)
        ω = @. (((1:J)) * piHigh / J)
        z = ω * 2^m
        phizJ, phizdJ = evaluateLogCharFuncAndDerivative(cf, -z[end], τ)
        phiz0, phizd0 = evaluateLogCharFuncAndDerivative(cf, 0.0, τ)
        phizdJ = -phizdJ
        phizd0 = -phizd0
        factor = 2^(m / 2) / J
        for k = 1:κ
            cl[k] -= factor / (12 * J) * real((1im * piHigh * k + 2^m * piHigh * phizdJ) * exp(1im * piHigh * k + phizJ) - (2^m * piHigh * phizd0) * exp(phiz0))
        end
        for k = -κ+1:0
            cl[k+2κ] -= factor / (12 * J) * real((1im * piHigh * k + 2^m * piHigh * phizdJ) * exp(1im * piHigh * k + phizJ) - (2^m * piHigh * phizd0) * exp(phiz0))
        end
    end
    return cl
end

struct VietaSlowDensityCalculator <: DensityCoeffCalculator
end

function calculate(calc::VietaSlowDensityCalculator, cf::CharFunc{MAINT,CR}, τ::T, m::Int, κ::Int, J::Int) where {T,MAINT,CR}
    piHigh = const_pi(cf)
    ω = @. (((1:J) - 1 / 2) * piHigh / J)
    z = ω * 2^m
    N = 2J
    phiz = zeros(CR, J)
    for (j, zj) = enumerate(z)
        phiz[j] = evaluateCharFunc(cf, -zj, τ) #equivalently w could evaluate at zj and multiply by exp(-ikw) instead of exp(ikw)
    end
    cl = zeros(T, 2 * κ)
    factor = 2^(m / 2) / J
    for k = 1:κ
        cl[k] = sum([real(phizj *exp(complex(zero(T), k * ωj))) for (phizj, ωj) = zip(phiz, ω)]) * factor
    end
    for k = -κ+1:0
        cl[k+2κ] = sum([phizj *real(exp(complex(zero(T), k * ωj))) for (phizj, ωj) = zip(phiz, ω)]) * factor
    end

    return cl
end

abstract type PayoffCoeffCalculator end

struct VietaPayoffCalculator <: PayoffCoeffCalculator
end

function calculate(calculator::VietaPayoffCalculator, p::SwiftCharFuncPricer{T}, strike, f, J, useVaryingJ) where {T}
    uk = p.uk
    m = p.m
    a = -p.c
    κ = p.k1
    ea = exp(a)
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    N = 2J
    ω = @. (((1:J) - 1 / 2) * piHigh * twom / J)
    payoffT = zeros(Complex{T}, N)
    #Gm,k = sum -eiwk payoffT[j] and we have Gm,1-k = e-ipi(k-1)/N fft(sum(eiwk payoffT[j]))
    #Gm,1+k = eipi(k+1)/N fft(sum(-eiwk payoffT[j]))
    #Fm,k = Gm,-k.   Fm,1-k = e-ipi(k-1)/N fft(sum(eiwk payoffT[j])), Gm,k-1 = e-ipi(k-1)/N fft(sum(eiwk payoffT[j]))
    @inbounds for j = 1:J
        payoffT[j] = (exp(ω[j] * logStrike * 1im) * (estrike / ω[j] - estrike * (ω[j] + 1im) / ((ω[j])^2 + one(T))) - exp(ω[j] * a * 1im) * (estrike / ω[j] - ea * (ω[j] + 1im) / ((ω[j])^2 + one(T)))) * (-1im)
    end
    fft!(payoffT)
    ukz = payoffT
    for k = 1:κ
        # s, c = sincos(- piHigh * (k - 1) / N); uk[k-1] = (real(ukz[k])* c -imag(ukz[k])*s) / J
        uk[k] = real(ukz[k+1] * exp(-1im * piHigh * k / N)) / J
    end
    for k = -κ+2:0
        uk[2κ+k-1] = real(ukz[N+k] * exp(-1im * piHigh * (k - 1) / N)) / J
    end
    uk[2κ] = real(ukz[1]) / J
end


struct VietaInfPayoffCalculator <: PayoffCoeffCalculator
end

function calculate(calculator::VietaInfPayoffCalculator, p::SwiftCharFuncPricer{T}, strike, f, J, useVaryingJ) where {T}
    uk = p.uk
    m = p.m
    κ = p.k1
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    #  a = -2p.c
    a = (-p.k2+1)/twom
    ea = exp(a)
    twom2 = 2^(m/2)
    N = 2J
    ω = @. (((1:J) - 1 / 2) * piHigh * twom / J)
    payoffT = zeros(Complex{T}, N)
    @inbounds for j = 1:J
        payoffT[j] = one(T)/(ω[j]*(one(T)+1im*ω[j]))*(exp((one(T)+1im*ω[j])*logStrike) -exp((one(T)+1im*ω[j])*a))
    end
    fft!(payoffT)
    ukz = payoffT
    for k = 1:κ
        uk[k] = imag(ukz[k+1] * exp(-1im * piHigh * k / N)) / J  + (estrike-ea)/(2*twom)
    end
    for k = -κ+2:0
        uk[2κ+k-1] = imag(ukz[N+k] * exp(-1im * piHigh * (k - 1) / N)) / J + (estrike-ea)/(2*twom)
    end
    uk[2κ] = imag(ukz[1]) / J + (estrike-ea)/(2*twom)
end

struct VietaDigitalPayoffCalculator <: PayoffCoeffCalculator
end

function calculate(calculator::VietaDigitalPayoffCalculator, p::SwiftCharFuncPricer{T}, strike, f, J) where {T}
    uk = p.uk
    m = p.m
    a = -p.c
    κ = p.k1
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    twom2 = 2^(m/2)
    N = 2J
    ω = @. (((1:J) - 1 / 2) * piHigh * twom / J)
    payoffT = zeros(Complex{T}, N)
    #Gm,k = sum -eiwk payoffT[j] and we have Gm,1-k = e-ipi(k-1)/N fft(sum(eiwk payoffT[j]))
    #Gm,1+k = eipi(k+1)/N fft(sum(-eiwk payoffT[j]))
    #Fm,k = Gm,-k.   Fm,1-k = e-ipi(k-1)/N fft(sum(eiwk payoffT[j])), Gm,k-1 = e-ipi(k-1)/N fft(sum(eiwk payoffT[j]))
    @inbounds for j = 1:J
        payoffT[j] = 1/ω[j]*exp(1im*ω[j]*logStrike)
    end
    fft!(payoffT)
    ukz = payoffT
    for k = 1:κ
        # s, c = sincos(- piHigh * (k - 1) / N); uk[k-1] = (real(ukz[k])* c -imag(ukz[k])*s) / J
        uk[k] = imag(ukz[k+1] * exp(-1im * piHigh * k / N)) / J  + 1/(2*twom)
    end
    for k = -κ+2:0
        uk[2κ+k-1] = imag(ukz[N+k] * exp(-1im * piHigh * (k - 1) / N)) / J + 1/(2*twom)
    end
    uk[2κ] = imag(ukz[1]) / J + 1/(2*twom)
end


struct EulerMaclaurinSlowPayoffCalculator{Order} <: PayoffCoeffCalculator #first kind of Euler-Maclaurin with or without derivative correction
end
struct SimpsonEulerMaclaurinPayoffCalculator{Order} <: PayoffCoeffCalculator #first kind of Euler-Maclaurin with or without derivative correction
end
#TODO yet another possiblity is to integrate teh digital put payoff from 0 to K. to obtain the vanilla put.

struct EulerMaclaurinPayoffCalculator{Order} <: PayoffCoeffCalculator 
end


function calculate(calculator::EulerMaclaurinPayoffCalculator{Order}, p::SwiftCharFuncPricer{T}, strike, f, J, useVaryingJ) where {T,Order}
    uk = p.uk
    m = p.m
    a = -p.c
    κ = p.k1
    ea = exp(a)
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    N = 2J
    ω = @. (((1:J)) * piHigh * twom / J)
    payoffT = zeros(Complex{T}, N)
    @inbounds for j = 1:J
        #payoffT[j+1] = -1im * (estrike*exp(1im* ω[j] * logStrike) - estrike*exp(1im* ω[j] * a))/ω[j] - (one(T) - 1im*ω[j])*(estrike*exp(1im* ω[j] * logStrike)-ea*exp(1im* ω[j] * a))/(ω[j]^2+one(T))
        payoffT[j+1] = (exp(ω[j] * logStrike * 1im) * (estrike / ω[j] - estrike * (ω[j] + 1im) / ((ω[j])^2 + one(T))) - exp(ω[j] * a * 1im) * (estrike / ω[j] - ea * (ω[j] + 1im) / ((ω[j])^2 + one(T)))) * (-1im)
    end
    payoffJ = payoffT[J+1]
    payoffT[1] = (estrike*(logStrike-a)-(estrike-ea))/2
    payoffT[J+1] /= 2 #necessary as trapezoidal is on parseval identity, not an expansion.
    fft!(payoffT)
    ukz = payoffT
    for k = 1:κ
        uk[k] = real(ukz[k+1])/J
    end
    for k = -κ+2:0
        uk[2κ+k-1] = real(ukz[N+k]) / J
    end
    uk[2κ] = real(ukz[1]) / J
    if Order == 1
        factor = piHigh/(12*J^2)
        term = (-1im * (logStrike*estrike*exp(1im* ω[J] * logStrike) - a*estrike*exp(1im* ω[J] * a))/ω[J] - (one(T) - 1im*ω[J])*(logStrike*estrike*exp(1im* ω[J] * logStrike)-a*ea*exp(1im* ω[J] * a))/(ω[J]^2+one(T)))
        term -= -estrike*(exp(1im* ω[J] * logStrike)-exp(1im* ω[J] * a))/ω[J]^2 - (one(T) - 1im*ω[J])^2* (estrike*exp(1im* ω[J] * logStrike)-ea*exp(1im* ω[J] * a))/(ω[J]^2+one(T))^2 
        term *= twom  
        for k = 1:κ
            uk[k] += imag(-k*payoffJ + term)*factor * (-one(T))^k
        end
        for k = -κ+1:0
            uk[2κ+k]+= imag(-k*payoffJ + term)*factor * (-one(T))^k
        end    
    end
end
function calculate(calculator::EulerMaclaurinSlowPayoffCalculator{N}, p::SwiftCharFuncPricer{T}, strike, f, J, useVaryingJ) where {T,N}
    uk = p.uk
    m = p.m
    a = -p.c
    κ = p.k1
    ea = exp(a)
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    for k = -κ+1:κ
        ki = if k > 0
            k
        else
            2κ + k
        end
        if useVaryingJ
            ln2J = ceil(Int, log2(twom * p.c + abs(k)))    # the paper advise log2(pi*κ) because of particular series expansion. log2(κ) is sufficient for FFT.    
            J = 2^ln2J * 2
        end
        ω = @. (((1:J)) * piHigh / J)
        weights = ones(T, length(ω))
        weights[end] /= 2
        trap = ((estrike * (logStrike - a) - (estrike - ea)) / 2 + sum([wj * real(exp(-ωj * k * 1im) * ((exp(ωj * twom * logStrike * 1im) - exp(ωj * twom * a * 1im)) * estrike / (ωj * twom * 1im) - (exp(ωj * twom * logStrike * 1im + logStrike) - exp(ωj * twom * a * 1im + a)) / (ωj * twom * 1im + one(T)))) for (wj, ωj) = zip(weights, ω)])) / J
        if N == 1
            integrand = function (y)
                (twom * y - k) * (estrike - exp(y)) * sin(piHigh * (twom * y - k))
            end
            trap += piHigh / (12 * J^2) * integrateSimpsonGG(integrand, a, logStrike, 1e-15, maxRecursionDepth=22)

            #todo rederive explicit formula using integration by parts and formula for int estrike-ey sin.
            # cim = piHigh * twom
            # sina, cosa = sincos(cim * a)
            # sinb, cosb = sincos(cim * logStrike)
            # is = -estrike*(cosb+cim*sinb)/(cim*(1+cim*cim)) + estrike*cosa/cim + ea*(sina-cim*cosa)/(1+cim*cim)
            # cim2 = cim^2
            # cim3 = cim^3
            # cim4 = cim^4
            # d = -estrike*(((cim4+cim2)*logStrike-3*cim2-1)*sinb+((cim3+cim)*logStrike+2*cim3)*cosb+(cim4+2*cim2+1)*sina-a*cim*(cim4+2*cim2+1)*cosa) - ea*((-(1+a)*cim4+(1-a)*cim2)*sina+(a*cim4+(a-2)*cim2)*cim*cosa)
            # d /= -(cim * (cim4 + cim2*2 + 1))
            # denom = (12 * J^2)
            # ratios = piHigh / denom * is
            # ratiod = d / denom
            #         trap -= (-1)^k * ( k  * ratios + ratiod)

        end
        uk[ki] = trap
    end
    #seems much worse than midpoint on Jp.
end

function calculate(calculator::SimpsonEulerMaclaurinPayoffCalculator{0}, p::SwiftCharFuncPricer{T}, strike, f, Jp, useVaryingJ) where {T}
    uk = p.uk
    m = p.m
    a = -p.c
    κ = p.k1
    ea = exp(a)
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    J = Jp >> 1
    for k = -κ+1:κ
        ki = if k > 0
            k
        else
            2κ + k
        end
        if useVaryingJ
            ln2J = ceil(Int, log2(abs(k)))    # the paper advise log2(pi*κ) because of particular series expansion. log2(κ) is sufficient for FFT.    
            J = 2^ln2J
        end
        ω = @. (((1:J)) * piHigh / J)
        weights = ones(T, length(ω))
        weights[end] /= 2
        trap = ((estrike * (logStrike - a) - (estrike - ea)) / 2 + sum([wj * real(exp(-ωj * k * 1im) * ((exp(ωj * twom * logStrike * 1im) - exp(ωj * twom * a * 1im)) * estrike / (ωj * twom * 1im) - (exp(ωj * twom * logStrike * 1im + logStrike) - exp(ωj * twom * a * 1im + a)) / (ωj * twom * 1im + one(T)))) for (wj, ωj) = zip(weights, ω)])) / J
        ω = @. (((1:J) - 1 / 2) * piHigh / J)
        midp = (sum([real(exp(-ωj * k * 1im) * ((exp(ωj * twom * logStrike * 1im) - exp(ωj * twom * a * 1im)) * estrike / (ωj * twom * 1im) - (exp(ωj * twom * logStrike * 1im + logStrike) - exp(ωj * twom * a * 1im + a)) / (ωj * twom * 1im + one(T)))) for ωj = ω])) / J
        uk[ki] = (midp + 2 * midp) / 3
    end
    #seems much worse than midpoint on Jp.
end

@enum Rule begin
    MidPoint
    Trapezoidal
    Simpson
    Boole
end
struct IntegralPayoffCalculator <: PayoffCoeffCalculator
    rule::Rule
end

function calculate(calculator::IntegralPayoffCalculator, p::SwiftCharFuncPricer{T}, strike, f, J, useVaryingJ) where {T}
    uk = p.uk
    m = p.m
    a = -p.c
    κ = p.k1
    ea = exp(a)
    piHigh = pi #  piHigh = const_pi(cf)
    logStrike = log(strike / f)
    estrike = strike / f
    twom = 2^m
    for k = -κ+1:κ
        integrand = function (y)
            (estrike - exp(y)) * sinc((twom * y - k))
        end
        ki = if k > 0
            k
        else
            2κ + k
        end
        if useVaryingJ
            ln2J = ceil(Int, log2(abs(k)))    # the paper advise log2(pi*κ) because of particular series expansion. log2(κ) is sufficient for FFT.    
            J = 2^ln2J
        end
        uk[ki] = if calculator.rule == MidPoint
            mid(integrand, a, logStrike, J)
        elseif calculator.rule == Simpson
            simpson(integrand, a, logStrike, J)
        elseif calculator.rule == Boole
            boole(integrand, a, logStrike, J)
        elseif calculator.rule == Trapezoidal
            trap(integrand, a, logStrike, J)
        end
    end
end

function priceEuropean(
    p::SwiftCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T; payoffCalculator=VietaPayoffCalculator(), J=p.j, useVaryingJ=false
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    local pricePut
    x = log(forward / strike)
    if x >= p.c
        pricePut = 0
    elseif x <= -p.c
        pricePut = discountDf * (strike - forward)
    else
        uk = p.uk
        m = p.m
        f = forward

        twom2 = 2^(m / 2)
        calculate(payoffCalculator, p, strike, f, J, useVaryingJ)
        sumPut = dot(p.cl, uk)
        pricePut = discountDf * f * sumPut * twom2
    end
    if isCall
        return pricePut + discountDf * (forward - strike)
    end
    return pricePut
end

function priceDigital(
    p::SwiftCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T;
    payoffCalculator=VietaDigitalPayoffCalculator(), J=p.j,
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    local pricePut
    x = log(forward / strike)
    if x >= p.c
        pricePut = 0
    elseif x <= -p.c
        pricePut = discountDf
    else   
        uk = p.uk
        m = p.m
        f = forward

        twom2 = 2^(m / 2)
        calculate(payoffCalculator, p, strike, f, J)
        sumPut = dot(p.cl, uk)
        pricePut = discountDf * sumPut * twom2
    end
    if isCall
        return -pricePut + discountDf
    end
    return pricePut
end


