export FCCCharFuncPricer

#Adaptive Filon-Clenshaw-Curtis quadrature of the Lewis formula.
#The amplitude g(u) = φ(u - i/2) / (u² + 1/4) is approximated on [0, b] by piecewise Chebyshev interpolants of degree N,
#on panels bisected until the last Chebyshev coefficients are negligible. The panels do not depend on the strike:
#for each strike, the oscillatory factor exp(-i u log(K/F)) is integrated exactly against the interpolants (Filon).
struct FCCCharFuncPricer{T,CR}
    τ::T
    b::T
    centers::Vector{T}
    halfWidths::Vector{T}
    coeffs::Matrix{CR} #Chebyshev coefficients, one column per panel
    vControl::T
    pi::T
    xGL::Vector{T} #Gauss-Legendre nodes, for panels where the oscillation is weak
    wvGL::Matrix{CR} #Gauss-Legendre weights times the interpolant values at the nodes, one column per panel
end

function FCCCharFuncPricer(
    p::CharFunc{MAINT,CR},
    τ::T;
    qTol::T=T(1e-10),
    b::T=zero(T),
    N::Int=16,
    maxDepth::Int=30,
) where {MAINT,CR,T}
    iPure = oneim(p)
    g = function (u)
        evaluateCharFunc(p, u - iPure / 2, τ) / (u^2 + one(T) / 4)
    end
    if b == zero(T)
        b = fccTruncation(g, qTol)
    end
    tNodes = [cospi(T(j) / N) for j = 0:N]
    centers = T[]
    halfWidths = T[]
    coeffList = Vector{Vector{typeof(g(zero(T)))}}()
    #depth-first bisection, panels are accepted from left to right
    stack = [(zero(T), b, 0)]
    while !isempty(stack)
        (a1, b1, depth) = pop!(stack)
        c = (a1 + b1) / 2
        h = (b1 - a1) / 2
        values = [g(c + h * t) for t in tNodes]
        coeff = chebCoefficientsFromValues(values)
        tailError = h * (abs(coeff[end]) + abs(coeff[end-1]))
        if tailError <= qTol * h / b || depth >= maxDepth
            push!(centers, c)
            push!(halfWidths, h)
            push!(coeffList, coeff)
        else
            push!(stack, (c, b1, depth + 1))
            push!(stack, (a1, c, depth + 1))
        end
    end
    coeffs = reduce(hcat, coeffList)
    xGL, wGL = gausslegendre(2N + 16)
    xGL = T.(xGL)
    wGL = T.(wGL)
    tGL = [cos(k * acos(x)) for x in xGL, k = 0:N]
    wvGL = wGL .* (tGL * coeffs)
    return FCCCharFuncPricer{T,eltype(coeffs)}(τ, b, centers, halfWidths, coeffs, T(getControlVariance(p)), const_pi(p), xGL, wvGL)
end

Base.broadcastable(p::FCCCharFuncPricer) = Ref(p)

#power of two b beyond which the integrand tail is negligible (checked on two consecutive points, as the integrand may oscillate).
function fccTruncation(g, tol::T) where {T}
    b = one(T)
    while b < T(1e8) && (abs(g(b)) * b > tol || abs(g(2b)) * 2b > tol)
        b *= 2
    end
    return 2b
end

#Chebyshev coefficients of the interpolant at the points cos(jπ/N), j=0..N.
function chebCoefficientsFromValues(values::AbstractVector)
    N = length(values) - 1
    coeff = similar(values)
    for k = 0:N
        s = (values[1] + values[N+1] * (-1)^k) / 2
        for j = 1:N-1
            s += values[j+1] * cospi(j * k / N)
        end
        coeff[k+1] = 2 * s / N
    end
    coeff[1] /= 2
    coeff[N+1] /= 2
    return coeff
end

#moments M_k = ∫_{-1}^{1} T_k(t) exp(iκt) dt, k=0..N, by forward recurrence, stable for N <= |κ|.
function chebOscillatoryMoments!(M::AbstractVector, κ::T) where {T}
    N = length(M) - 1
    #forward recurrence, stable for k <= |κ|. Uses ∫T'_m e^{iκt} = B_m - iκ M_m with B_m = e^{iκ} - (-1)^m e^{-iκ}
    #and 2 T_k = T'_{k+1}/(k+1) - T'_{k-1}/(k-1).
    eip = cis(κ)
    eim = conj(eip)
    bnd(m) = isodd(m) ? eip + eim : eip - eim
    iκ = im * κ
    M[1] = 2 * sin(κ) / κ
    N >= 1 && (M[2] = (bnd(1) - M[1]) / iκ)
    N >= 2 && (M[3] = (bnd(2) - 4 * M[2]) / iκ)
    for k = 2:N-1
        M[k+2] = (bnd(k + 1) - (k + 1) * (2 * M[k+1] + (bnd(k - 1) - iκ * M[k]) / (k - 1))) / iκ
    end
    return M
end

function priceEuropean(
    p::FCCCharFuncPricer{T,CR},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T,CR}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
    end
    ω = log(strike / forward)
    N = size(p.coeffs, 1) - 1
    M = zeros(Complex{T}, N + 1)
    integral = zero(T)
    @inbounds for i = eachindex(p.centers)
        h = p.halfWidths[i]
        κ = -ω * h
        s = zero(Complex{T})
        if abs(κ) >= N
            chebOscillatoryMoments!(M, κ)
            for k = 1:N+1
                s += p.coeffs[k, i] * M[k]
            end
        else
            #weak oscillation: Gauss-Legendre on the interpolant
            for j = eachindex(p.xGL)
                s += p.wvGL[j, i] * cis(κ * p.xGL[j])
            end
        end
        integral += real(h * cis(-ω * p.centers[i]) * s)
    end
    callPrice = -discountDf * sqrt(strike * forward) / p.pi * integral
    if p.vControl != 0
        #use control variate
        sqrtVar = sqrt(p.vControl * p.τ)
        d1 = log(forward / strike) / sqrtVar + sqrtVar / 2
        d2 = d1 - sqrtVar
        callPrice += discountDf * (forward * normcdf(d1) - strike * normcdf(d2))
    else
        callPrice += discountDf * forward
    end
    if isCall
        return callPrice
    end
    return callPrice - (forward - strike) * discountDf
end
