export FCCCharFuncPricer

#Adaptive Filon-Clenshaw-Curtis quadrature of the Lewis formula.
#The amplitude g(u) = φ(u - i/2) / (u² + 1/4) is approximated on [0, b] by piecewise Chebyshev interpolants,
#on panels refined until the last Chebyshev coefficients are negligible. The panels do not depend on the strike:
#for each strike, the oscillatory factor exp(-i u log(K/F)) is integrated exactly against the interpolants (Filon).
#With hp=true, a panel which is not resolved at degree N is first tried at degree 2N (the Clenshaw-Curtis points are
#nested, only N new evaluations are needed) before being bisected.
struct FCCCharFuncPricer{T,CR}
    τ::T
    b::T
    centers::Vector{T}
    halfWidths::Vector{T}
    coeffs::Vector{Vector{CR}} #Chebyshev coefficients of each panel
    level::Vector{Int} #1 for degree N, 2 for degree 2N
    vControl::T
    pi::T
    xGL::Vector{Vector{T}} #Gauss-Legendre nodes of each level, for panels where the oscillation is weak
    wvGL::Vector{Vector{CR}} #Gauss-Legendre weights times the interpolant values at the nodes, for each panel
end

function FCCCharFuncPricer(
    p::CharFunc{MAINT,CR},
    τ::T;
    qTol::T=T(1e-10),
    b::T=zero(T),
    N::Int=16,
    hp::Bool=true,
    dyadic::Bool=true, #initial panels [0,1], [1,2], [2,4], ..., [b/2,b], otherwise nInitial equal panels
    nInitial::Int=4,
    maxDepth::Int=30,
    maxEvaluations::Int=100_000,
) where {MAINT,CR,T}
    iPure = oneim(p)
    gRaw = function (u)
        evaluateCharFunc(p, u - iPure / 2, τ) / (u^2 + one(T) / 4)
    end
    VT = typeof(gRaw(zero(T)))
    #panel endpoints and midpoints are shared between panels (and with the truncation search)
    cache = Dict{T,VT}()
    g(u) = get!(() -> gRaw(u), cache, u)
    if b == zero(T)
        b = fccTruncation(g, qTol)
    end
    rule = fccRule(T, N)
    tNodes, tNodes2, C1, C2 = rule.tNodes, rule.tNodes2, rule.C1, rule.C2
    #Clenshaw-Curtis nodes of [a1, b1], the end and mid points exactly
    node(a1, b1, c, h, t, j, n) = j == 0 ? b1 : (j == n ? a1 : (2j == n ? c : c + h * t))
    centers = T[]
    halfWidths = T[]
    coeffList = Vector{Vector{VT}}()
    levels = Int[]
    #resolved when the Chebyshev tail is below the tolerance or at the level of the round-off in the values
    isResolved(coeff, values, h) = begin
        tail = abs(coeff[end]) + abs(coeff[end-1])
        tail <= qTol / b || tail <= 64 * length(coeff) * eps(T) * maximum(abs, values)
    end
    #depth-first bisection, panels are accepted from left to right
    bp = if dyadic
        vcat(zero(T), [T(2)^k for k = 0:floor(Int, log2(b)) if T(2)^k < b], b)
    else
        collect(range(zero(T), b, length=nInitial + 1))
    end
    stack = [(bp[i], bp[i+1], 0) for i = length(bp)-1:-1:1]
    while !isempty(stack)
        (a1, b1, depth) = pop!(stack)
        c = (a1 + b1) / 2
        h = (b1 - a1) / 2
        values = [g(node(a1, b1, c, h, tNodes[j+1], j, N)) for j = 0:N]
        if !all(isfinite, values)
            throw(DomainError(b1, string("the characteristic function is not finite on [", a1, ", ", b1, "], use a smaller truncation b")))
        end
        coeff = C1 * values
        level = 1
        resolved = isResolved(coeff, values, h)
        if !resolved && hp
            values2 = Vector{VT}(undef, 2N + 1)
            values2[1:2:end] = values
            for j = 2:2:2N
                values2[j] = g(node(a1, b1, c, h, tNodes2[j], j - 1, 2N))
            end
            if !all(isfinite, values2)
                throw(DomainError(b1, string("the characteristic function is not finite on [", a1, ", ", b1, "], use a smaller truncation b")))
            end
            coeff = C2 * values2
            level = 2
            resolved = isResolved(coeff, values2, h)
        end
        if !resolved && length(cache) >= maxEvaluations
            @warn "FCC: maximum number of evaluations reached, the tolerance may not be met" maxEvaluations
            resolved = true
        end
        if resolved || depth >= maxDepth
            push!(centers, c)
            push!(halfWidths, h)
            push!(coeffList, coeff)
            push!(levels, level)
        else
            push!(stack, (c, b1, depth + 1))
            push!(stack, (a1, c, depth + 1))
        end
    end
    wvGL = [rule.wtGL[l] * coeff for (coeff, l) in zip(coeffList, levels)]
    return FCCCharFuncPricer{T,VT}(τ, b, centers, halfWidths, coeffList, levels, T(getControlVariance(p)), const_pi(p), rule.xGL, wvGL)
end

Base.broadcastable(p::FCCCharFuncPricer) = Ref(p)

#Clenshaw-Curtis points and coefficient matrices for degrees N and 2N, Gauss-Legendre nodes and (weights times Chebyshev
#polynomials at the nodes) for each degree. They do not depend on the model and are computed once per (T, N).
const FCC_RULES = Dict{Tuple{DataType,Int},Any}()
const FCC_RULES_LOCK = ReentrantLock()
function fccRule(::Type{T}, N::Int) where {T}
    lock(FCC_RULES_LOCK) do
        get!(FCC_RULES, (T, N)) do
            xGL = Vector{Vector{T}}(undef, 2)
            wtGL = Vector{Matrix{T}}(undef, 2)
            for (l, n) in enumerate((N, 2N))
                x, w = gausslegendre(2n + 16)
                xGL[l] = T.(x)
                wtGL[l] = T.(w) .* [cos(k * acos(xi)) for xi in xGL[l], k = 0:n]
            end
            (tNodes=[cospi(T(j) / N) for j = 0:N], tNodes2=[cospi(T(j) / (2N)) for j = 0:2N],
                C1=chebCoefficientMatrix(T, N), C2=chebCoefficientMatrix(T, 2N), xGL=xGL, wtGL=wtGL)
        end
    end::NamedTuple{(:tNodes, :tNodes2, :C1, :C2, :xGL, :wtGL),Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T},Vector{Vector{T}},Vector{Matrix{T}}}}
end

#power of two b beyond which the integrand tail is negligible (checked on two consecutive points, as the integrand may oscillate).
#Stops at the last point where the integrand is finite (e.g. explicit numerical schemes for the characteristic function).
function fccTruncation(g, tol::T) where {T}
    b = one(T)
    while b < T(1e8)
        g1, g2 = g(b), g(2b)
        if !isfinite(g2)
            @warn "FCC: the characteristic function is not finite beyond u = $b, truncating there"
            return b
        end
        if abs(g1) * b <= tol && abs(g2) * 2b <= tol
            break
        end
        b *= 2
    end
    return 2b
end

#matrix mapping the values at the points cos(jπ/N), j=0..N, to the Chebyshev coefficients of the interpolant (DCT-I).
function chebCoefficientMatrix(::Type{T}, N::Int) where {T}
    C = zeros(T, N + 1, N + 1)
    for k = 0:N, j = 0:N
        w = (j == 0 || j == N) ? one(T) / 2 : one(T)
        C[k+1, j+1] = 2 * w * cospi(T(j * k) / N) / N
    end
    C[1, :] ./= 2
    C[N+1, :] ./= 2
    return C
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
    M = zeros(Complex{T}, maximum(length, p.coeffs))
    integral = zero(T)
    @inbounds for i = eachindex(p.centers)
        h = p.halfWidths[i]
        κ = -ω * h
        coeff = p.coeffs[i]
        n = length(coeff) - 1
        s = zero(Complex{T})
        if abs(κ) >= n
            Mi = view(M, 1:n+1)
            chebOscillatoryMoments!(Mi, κ)
            for k = 1:n+1
                s += coeff[k] * Mi[k]
            end
        else
            #weak oscillation: Gauss-Legendre on the interpolant
            x = p.xGL[p.level[i]]
            wv = p.wvGL[i]
            for j = eachindex(x)
                s += wv[j] * cis(κ * x[j])
            end
        end
        integral += real(h * cis(-ω * p.centers[i]) * s)
    end
    return lewisPrice(isCall, integral, strike, forward, p.τ, discountDf, p.vControl, p.pi)
end
