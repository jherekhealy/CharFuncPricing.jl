import DataStructures: SortedDict
using SpecialFunctions

normcdf(z::T) where {T} = erfc(-z / sqrt(T(2))) / 2

export FlinnCharFuncPricer, priceEuropean, AdaptiveFlinnCharFuncPricer


getControlVariance(
    ::Union{CharFunc{HestonParams{T}},CharFunc{SchobelZhuParams{T}}},
    τ::T,
) where {T} = T(0)

@inline function getControlVariance(
    cf::CVCharFunc{HestonParams{T},BlackParams{T}},
    τ::T,
)::T where {T}
    p = model(cf)
    # ektm = 1-exp(-p.κ*τ)
    # return ektm*(p.v0-p.θ)/(p.κ*τ) + p.θ
    p.v0
end


#Adaptive Flinn with variable transform.
struct AdaptiveFlinnCharFuncPricer{T}
    τ::T
    kcos::Array{T,2}
    vControl::T
    pi::T
    function AdaptiveFlinnCharFuncPricer(
        p::CharFunc,
        τ::T;
        qTol::T = sqrt(eps(T)),
    ) where {T}
        mcos = Dict{T,Tuple{T,T,T,T}}()
        iPure = oneim(p)
        @inline function integrand(z::T)::Tuple{T,T,T,T} where {T}
            if z == one(T)
                #x=+infty
                return zero(T), zero(T), zero(T), zero(T)
            end
            x = (1 + z) / (1 - z)
            u = x - iPure / 2
            (phi, phi_d) = evaluateCharFuncAndDerivative(p, u, τ) #derivative towards real part of u
            denom = (1 / 4 + x^2)
            phi /= denom
            phi_d = (phi_d - phi * 2 * x) / denom
            factor = 2 / (1 - z)^2
            phi *= factor
            phi_d = phi_d * factor^2 + 2 / (1 - z) * phi
            return real(phi), imag(phi), real(phi_d), imag(phi_d)
        end

        @inline function fcos(x::T)::Tuple{T,T} where {T}
            v = integrand(x)
            if (x != 1)
                mcos[x] = v
            end
            return (v[1], v[3])
        end
        @inline function fsin(x::T)::Tuple{T,T} where {T}
            # v = get!(mcos,x , integrand(x))
            local v
            if haskey(mcos, x)
                v = mcos[x]
            else
                v = integrand(x)
                if (x != 1)
                    mcos[x] = v
                end
            end
            return (v[2], v[4])
        end
        integrateQuinticHermite(fcos, -one(T), one(T), qTol, 16)
        integrateQuinticHermite(fsin, -one(T), one(T), qTol, 16)
        sortedKeys = sort!(collect(keys(mcos)))
        kcos = zeros(T, (5, length(mcos) * 2 - 1))
        #convert to original coordinate
        for (i, u) in enumerate(sortedKeys)
            v = mcos[u]
            x = (1 + u) / (1 - u)
            factor = 2 / (1 - u)^2
            kcos[:, i*2-1] = [
                x,
                v[1] / factor,
                v[2] / factor,
                (v[3] - 2 * v[1] / (1 - u)) / factor^2,
                (v[4] - 2 * v[2] / (1 - u)) / factor^2,
            ]
        end
        #add mid points
        for i = 1:length(mcos)-1
            a = @view kcos[:, i*2-1]
            b = @view kcos[:, (i+1)*2-1]
            x = (a[1] + b[1]) / 2
            u = (x - 1) / (x + 1) #(1-u)*x = 1+u => u(1+x) = x-1 => u=(x-1)/(x+1)
            v = integrand(u)
            factor = 2 / (1 - u)^2
            kcos[:, i*2] = [
                x,
                v[1] / factor,
                v[2] / factor,
                (v[3] - 2 * v[1] / (1 - u)) / factor^2,
                (v[4] - 2 * v[2] / (1 - u)) / factor^2,
            ]
        end
        #end interval = kcos[1,end]
        return new{T}(τ, kcos, getControlVariance(p, τ), const_pi(p))
    end

end


function priceEuropean(
    p::AdaptiveFlinnCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ",p.τ)))
    end
    #intervals are not equals anymore! => cubic hermite spline f, fp.
    #for now, double the points, note: need to transform f and f' back as well.
    freq = log(strike / forward) #-k
    sqrtEps = sqrt(eps(T))
    integral = zero(T)
    diff = zero(T)
    diff1 = zero(T)
    kcos = p.kcos
    nx = size(kcos)[2] #cols

    @inbounds for i = 1:2:nx-2
        fp = FlinnParams(freq, kcos[1, i], kcos[1, i+2])
        integral2 = flinnCosSin3Points(
            fp,
            kcos[1, i],
            kcos[1, i+1],
            kcos[1, i+2],
            kcos[2, i],
            kcos[2, i+1],
            kcos[2, i+2],
            kcos[3, i],
            kcos[3, i+1],
            kcos[3, i+2],
            kcos[4, i],
            kcos[4, i+1],
            kcos[4, i+2],
            kcos[5, i],
            kcos[5, i+1],
            kcos[5, i+2],
            freq,
        )
        integral += integral2
    end
    callPrice = -discountDf * sqrt(strike * forward) / p.pi * integral
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

#Flinn with truncation: very fast, but truncation rule not always robust
struct FlinnCharFuncPricer{T}
    τ::T
    b::T
    kcos::Array{T,2}
    vControl::T
    pi::T
    function FlinnCharFuncPricer(
        p::CharFunc,
        τ::T;
        tTol::T = T(1e-8),
        qTol::T = T(1e-8),
        b::T = Base.zero(T),
    ) where {T}
        if b == 0
            b = computeTruncation(p, τ, tTol)
        end
        mcos = Dict{T,Tuple{T,T,T,T}}()
        iPure = oneim(p)
        @inline function integrand(x::T)::Tuple{T,T,T,T} where {T}
            u = x - 0.5 * iPure
            (phi, phi_d) = evaluateCharFuncAndDerivative(p, u, τ) #derivative towards real part of u
            denom = 0.25 + x^2
            phi /= denom
            phi_d /= denom
            return real(phi),
            imag(phi),
            (real(phi_d) - real(phi) * 2 * x / denom),
            (imag(phi_d) - imag(phi) * 2 * x / denom)
        end
        @inline function fcos(x::T)::Tuple{T,T} where {T}
            v = integrand(x)
            mcos[x] = v
            return (v[1], v[3])
        end
        @inline function fsin(x::T)::Tuple{T,T} where {T}
            # v = get!(mcos,x , integrand(x))
            local v
            if haskey(mcos, x)
                v = mcos[x]
            else
                v = integrand(x)
                mcos[x] = v
            end
            return (v[2], v[4])
        end
        #TODO manual stack management.
        integrateQuinticHermite(fcos, T(0), b, qTol, 16)
        integrateQuinticHermite(fsin, T(0), b, qTol, 16)
        sortedKeys = sort!(collect(keys(mcos)))
        kcos = zeros(T, (5, length(mcos)))
        for (i, u) in enumerate(sortedKeys)
            v = mcos[u]
            kcos[:, i] = [u, v[1], v[2], v[3], v[4]]
        end
        return new{T}(τ, b, kcos, getControlVariance(p, τ), const_pi(p))
    end
end

function priceEuropean(
    p::FlinnCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    freq = log(strike / forward) #-k
    sqrtEps = sqrt(eps(T))
    integral = T(0)
    diff = T(0)
    diff1 = T(0)
    kcos = p.kcos
    nx = size(kcos)[2] #cols
    local fp::FlinnParams

    @inbounds for i = 1:2:nx-2
        if abs(diff - (kcos[1, i+2] - kcos[1, i])) > 1e-8
            fp = FlinnParams(freq, kcos[1, i], kcos[1, i+2])
            diff = kcos[1, i+2] - kcos[1, i]
        end
        integral2 = flinnCosSin3Points(
            fp,
            kcos[1, i],
            kcos[1, i+1],
            kcos[1, i+2],
            kcos[2, i],
            kcos[2, i+1],
            kcos[2, i+2],
            kcos[3, i],
            kcos[3, i+1],
            kcos[3, i+2],
            kcos[4, i],
            kcos[4, i+1],
            kcos[4, i+2],
            kcos[5, i],
            kcos[5, i+1],
            kcos[5, i+2],
            freq,
        )
        integral += integral2
    end
    callPrice = -discountDf * sqrt(strike * forward) / p.pi * integral
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


struct FlinnParams{T}
    M::T
    N::T
    P::T
    Q::T
    R::T
    S::T
    cost::T
    sint::T
end

function FlinnParams(t::T, x0::T, x2::T) where {T}
    theta = t * (x2 - x0) / 2
    return computeMNPQRS(theta)
end

function belowFlinnThreshold(theta::T)::Bool where {T}
    return abs(theta) <= 29*eps(T)^0.1 #0.8 for Float64
end

function belowFlinnThreshold(theta::Float64)::Bool
    return theta <= 0.8
end

@inline function computeMNPQRS(theta::T) where {T}
    theta2 = theta^2
    local M, N, P, Q, R, S
    sint, cost = sincos(theta)

    if belowFlinnThreshold(theta)
        M =
            theta * (
                -16.0 / 105.0 +
                theta2 * (
                    8.0 / 945 +
                    theta2 *
                    (-2.0 / 10395 + theta2 * (1.0 / 405405.0 - theta2 / 48648600.0))
                )
            )
        N =
            16.0 / 15.0 +
            theta2 * (
                -8.0 / 105 +
                theta2 * (2.0 / 945 + theta2 * (-1 / 31185.0 + theta2 / 3243240.0))
            )
        P =
            -1.0 / 15.0 +
            theta2 * (
                2.0 / 105 +
                theta2 *
                (-1.0 / 315 + theta2 * (2.0 / 7425 + theta2 * (-62.0 / 4729725.0)))
            )
        Q =
            theta * (
                -8 / 105.0 +
                theta2 * (
                    16.0 / 945 +
                    theta2 *
                    (-104.0 / 51975 + theta2 * (256.0 / 2027025 - theta2 * 16.0 / 3274425))
                )
            )
        R =
            14.0 / 15.0 +
            theta2 * (
                -16.0 / 105.0 +
                theta2 *
                (22.0 / 945 + theta2 * (-304.0 / 155925 + theta2 * (268.0 / 2837835)))
            )
        S =
            theta * (
                19.0 / 105 +
                theta2 * (
                    -2.0 / 63 +
                    theta2 * (1.0 / 275 + theta2 * (-2.0 / 8775 + theta2 * 34.0 / 3869775))
                )
            )
    else
        theta4 = theta2 * theta2
        theta6 = theta2 * theta4
        M = (16 * theta * (15 - theta2) * cost + 48 * (2 * theta2 - 5) * sint) / theta6
        N = (16 * theta * (3 - theta2) * sint - 48 * theta2 * cost) / theta6
        P =
            (
                2 * theta * (theta2 - 24) * sint * cost +
                15 * (theta2 - 4) * cost * cost +
                theta4 - 27 * theta2 + 60
            ) / theta6
        Q =
            (
                2 * (
                    theta * (12 - 5 * theta2) +
                    15 * (theta2 - 4) * sint * cost +
                    2 * theta * (24 - theta2) * cost * cost
                )
            ) / theta6
        R =
            (
                2 * (
                    theta * (156 - 7 * theta2) * sint * cost +
                    3 * (60 - 17 * theta2) * cost * cost - 15 * (12 - 5 * theta2)
                )
            ) / theta6
        S =
            (
                theta * (theta4 + 8 * theta2 - 24) +
                theta * (7 * theta2 - 156) * cost * cost +
                3 * (60 - 17 * theta2) * sint * cost
            ) / theta6
    end
    return FlinnParams(M, N, P, Q, R, S, cost, sint)
end


@inline function flinnCosSin3Points(
    p::FlinnParams{T},
    x0::T,
    x1::T,
    x2::T,
    fcos0::T,
    fcos1::T,
    fcos2::T,
    fsin0::T,
    fsin1::T,
    fsin2::T,
    fpcos0::T,
    fpcos1::T,
    fpcos2::T,
    fpsin0::T,
    fpsin1::T,
    fpsin2::T,
    t::T,
)::T where {T}
    if x0 == x2
        return T(0)
    end

    h = (x2 - x0) / 2

    st0, ct0 = sincos(t * x0)
    st1 = p.sint * ct0 + p.cost * st0
    ct1 = p.cost * ct0 - p.sint * st0
    st2 = p.sint * ct1 + p.cost * st1
    ct2 = p.cost * ct1 - p.sint * st1

    s2n = fsin0 * st0 / 2
    s2n += fsin2 * st2 / 2

    sp2n = fpsin0 * ct0 / 2
    sp2n += fpsin2 * ct2 / 2

    s2nm1 = fsin1 * st1
    sp2nm1 = fpsin1 * ct1

    value =
        h * (
            p.S * (-fsin2 * ct2 + fsin0 * ct0) +
            h * p.P * (fpsin2 * st2 - fpsin0 * st0) +
            p.R * s2n - h * p.Q * sp2n + p.N * s2nm1 - h * p.M * sp2nm1
        )
    c2n = fcos0 * ct0 / 2
    c2n += fcos2 * ct2 / 2

    cp2n = fpcos0 * st0 / 2
    cp2n += fpcos2 * st2 / 2

    c2nm1 = fcos1 * ct1
    cp2nm1 = fpcos1 * st1

    value +=
        h * (
            p.S * (fcos2 * st2 - fcos0 * st0) +
            h * p.P * (fpcos2 * ct2 - fpcos0 * ct0) +
            p.R * c2n +
            h * p.Q * cp2n +
            p.N * c2nm1 +
            h * p.M * cp2nm1
        )
    return value
end

computeTruncation(p::CVCharFunc, τ::Float64, tol::Float64) =
    computeTruncation(p.main, τ, tol)

function computeTruncation(cf::CharFunc{HestonParams{T}}, τ::T, tol::T) where {T}
    p = model(cf)
    c_inf = cinf(p, τ)
    u = T(lambertW(Float64(c_inf / tol))) / c_inf
    if p.v0 * τ < 0.1
        ushort = sqrt(T(lambertW(Float64(p.v0 * τ / (tol^2)))) / (p.v0 * τ))
        u = max(u, ushort)
    end
    #u = max(u,10.0)
    #end
    # println("umax ",u)
    return u
end

const p2 = (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
const h0 = (
    0.0762880721720079,
    -0.0329193509714785,
    -0.3438600114643695,
    0.5060431337241802,
    1.5888963130793203,
    0.5060431337241802,
    -0.3438600114643695,
    -0.0329193509714785,
    0.0762880721720079,
)
const h1 = (
    1.60554439048059e-03,
    -3.96727024402095e-02,
    -1.57110962304862e-01,
    -2.32981551813796e-01,
    0.0,
    2.32981551813796e-01,
    1.57110962304862e-01,
    3.96727024402095e-02,
    -1.60554439048059e-03,
)

function integrateQuinticHermite(
    f,
    a::T,
    b::T,
    tol::T,
    maxRecursionDepth::Int,
)::T where {T}
    if tol < eps(T)
        tol = eps(T)
    end

    if b <= a
        return Base.zero(T)
    end
    c = (a + b) / 2
    h = (b - a) / 2

    is = Base.zero(T)
    x = Vector{T}(undef, length(p2))
    y = Vector{T}(undef, length(p2))
    yp = Vector{T}(undef, length(p2))
    @inbounds for i = 1:length(p2)
        x[i] = c + h * p2[i]
        y[i], yp[i] = f(x[i])
        is += h * (h0[i] * y[i] + h * h1[i] * yp[i])
    end
    if is == Base.zero(T)
        is = one(T)
    end
    is = is * tol / eps(T)

    @inbounds result =
        adaptiveQuinticHermiteAux(
            f,
            x[1],
            x[3],
            is,
            y[1],
            y[3],
            y[2],
            yp[1],
            yp[3],
            yp[2],
            maxRecursionDepth,
        ) +
        adaptiveQuinticHermiteAux(
            f,
            x[3],
            x[5],
            is,
            y[3],
            y[5],
            y[4],
            yp[3],
            yp[5],
            yp[4],
            maxRecursionDepth,
        ) +
        adaptiveQuinticHermiteAux(
            f,
            x[5],
            x[7],
            is,
            y[5],
            y[7],
            y[6],
            yp[5],
            yp[7],
            yp[6],
            maxRecursionDepth,
        ) +
        adaptiveQuinticHermiteAux(
            f,
            x[7],
            x[9],
            is,
            y[7],
            y[9],
            y[8],
            yp[7],
            yp[9],
            yp[8],
            maxRecursionDepth,
        )
    return result
end

function adaptiveQuinticHermiteAux(
    f,
    a::T,
    b::T,
    is::T,
    fa::T,
    fb::T,
    fc::T,
    fpa::T,
    fpb::T,
    fpc::T,
    bottom::Int,
)::T where {T}
    c = (a + b) / 2
    h = b - a
    d = (a + c) / 2
    e = (c + b) / 2
    fd, fpd = f(d)
    fe, fpe = f(e)
    i1 = (h / 30) * (7 * fa + 16 * fc + 7 * fb + h / 2 * (fpa - fpb))
    Sleft = (h / 60) * (7 * fa + 16 * fd + 7 * fc + h / 4 * (fpa - fpc))
    Sright = (h / 60) * (7 * fc + 16 * fe + 7 * fb + h / 4 * (fpc - fpb))
    i2 = Sleft + Sright
    if ((is + (i1 - i2)/63) == is) || (d == e)
        return i2
    end
    if bottom <= 0
        # println("Number of maximum recursions reached")
        return i2
    end
    result =
        adaptiveQuinticHermiteAux(f, a, c, is, fa, fc, fd, fpa, fpc, fpd, bottom - 1) +
        adaptiveQuinticHermiteAux(f, c, b, is, fc, fb, fe, fpc, fpb, fpe, bottom - 1)
    return result
end
