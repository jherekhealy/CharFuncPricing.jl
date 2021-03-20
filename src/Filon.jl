import DataStructures: SortedDict

export priceEuropean, AdaptiveFilonCharFuncPricer


#Adaptive Filon with variable transform.
struct AdaptiveFilonCharFuncPricer{T}
    τ::T
    kcos::Array{T,2}
    vControl::T
    pi::T
    function AdaptiveFilonCharFuncPricer(
        p::CharFunc,
        τ::T;
        qTol::T = sqrt(eps(T)),
    ) where {T}
        mcos = Dict{T,Tuple{T,T}}()
        iPure = oneim(p)
        @inline function integrand(z::T)::Tuple{T,T} where {T}
            if z == one(T)
                #x=+infty
                return zero(T), zero(T)
            end
            x = (1 + z) / (1 - z)
            u = x - iPure / 2
            phi = evaluateCharFunc(p, u, τ)
            denom = (1 / 4 + x^2)
            phi /= denom
            factor = 2 / (1 - z)^2
            phi *= factor
            return real(phi), imag(phi)
        end

        @inline function fcos(x::T)::T where {T}
            v = integrand(x)
            if (x != 1)
                mcos[x] = v
            end
            return v[1]
        end
        @inline function fsin(x::T)::T where {T}
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
            return v[2]
        end
        integrateSimpsonGG(fcos, -one(T), one(T), qTol, 16)
        integrateSimpsonGG(fsin, -one(T), one(T), qTol, 16)
        sortedKeys = sort!(collect(keys(mcos)))
        kcos = zeros(T, (3, length(mcos) * 2 - 1))
        #convert to original coordinate
        for (i, u) in enumerate(sortedKeys)
            v = mcos[u]
            x = (1 + u) / (1 - u)
            factor = 2 / (1 - u)^2
            kcos[:, i*2-1] = [x, v[1] / factor, v[2] / factor]
        end
        #add mid points
        for i = 1:length(mcos)-1
            a = @view kcos[:, i*2-1]
            b = @view kcos[:, (i+1)*2-1]
            x = (a[1] + b[1]) / 2
            u = (x - 1) / (x + 1) #(1-u)*x = 1+u => u(1+x) = x-1 => u=(x-1)/(x+1)
            v = integrand(u)
            factor = 2 / (1 - u)^2
            kcos[:, i*2] = [x, v[1] / factor, v[2] / factor]
        end
        #end interval = kcos[1,end]
        return new{T}(τ, kcos, getControlVariance(p, τ), const_pi(p))
    end

end


function priceEuropean(
    p::AdaptiveFilonCharFuncPricer{T},
    isCall::Bool,
    strike::T,
    forward::T,
    τ::T,
    discountDf::T,
) where {T}
    if τ != p.τ
        throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
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
        fp = FilonParams(freq, kcos[1, i], kcos[1, i+2])
        integral2 = filonCosSin3Points(
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


struct FilonParams{T}
    α::T
    β::T
    γ::T
    cost::T
    sint::T
end

function FilonParams(t::T, x0::T, x2::T) where {T}
    theta = t * (x2 - x0) / 2
    return computeFilonParams(theta)
end

function belowFilonThreshold(theta::T)::Bool where {T}
    return abs(theta) <= 15 * eps(T)^0.125
end

function belowFilonThreshold(theta::Float64)::Bool
    return theta <= 1.0 / 6.0
end

@inline function computeFilonParams(theta::T) where {T}
    theta2 = theta^2
    local α, β, γ
    sint, cost = sincos(theta)

    if belowFilonThreshold(theta)
        α = theta * (theta2 * (2 / 45 + theta2 * (-2 / 315 + 2 * theta2 / 4725)))
        β =
            2 / 3 +
            theta2 *
            (2 / 15 + (theta2 * (-4 / 105 + theta2 * (2 / 567 - 4 * theta2 / 22275))))
        γ = 4 / 3 + theta2 * (-2 / 15 + theta2 * (1 / 210 - theta2 / 11340))
    else
        theta3 = theta2 * theta
        α = (theta2 + theta * sint * cost - 2.0 * sint * sint) / theta3
        β = (2 * theta + 2 * theta * cost * cost - 4 * sint * cost) / theta3
        γ = 4 * (sint - theta * cost) / theta3
    end
    return FilonParams(α, β, γ, cost, sint)
end


@inline function filonCosSin3Points(
    p::FilonParams{T},
    x0::T,
    x1::T,
    x2::T,
    fcos0::T,
    fcos1::T,
    fcos2::T,
    fsin0::T,
    fsin1::T,
    fsin2::T,
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

    s2nm1 = fsin1 * st1

    value = h * (p.α * (fsin0 * ct0 - fsin2 * ct2) + p.β * s2n + p.γ * s2nm1)
    c2n = fcos0 * ct0 / 2
    c2n += fcos2 * ct2 / 2

    c2nm1 = fcos1 * ct1
    value += h * (p.α * (fcos2 * st2 - fcos0 * st0) + p.β * c2n + p.γ * c2nm1)
    return value
end


function integrateSimpsonGG(f, a::T, b::T, tol::T, maxRecursionDepth::Int)::T where {T}
    if tol < eps(T)
        tol = eps(T)
    end

    if b <= a
        return Base.zero(T)
    end
    c = (a + b) / 2

    fa = f(a)
    fm = f(c)
    fb = f(b)
    yy0 = f(a + 0.9501 * (b - a))
    yy1 = f(a + 0.2311 * (b - a))
    yy2 = f(a + 0.6068 * (b - a))
    yy3 = f(a + 0.4860 * (b - a))
    yy4 = f(a + 0.8913 * (b - a))

    is = (b - a) / 8 * (fa + fm + fb + yy0 + yy1 + yy2 + yy3 + yy4)
    if is == 0
        is = b - a
    end
    is = is * tol / eps(T)

    result = adaptiveSimpsonAux(f, a, b, fa, fm, fb, is, maxRecursionDepth)
    return result
end

function adaptiveSimpsonAux(
    f,
    a::T,
    b::T,
    fa::T,
    fm::T,
    fb::T,
    is::T,
    bottom::Int,
)::T where {T}
    m = (b + a) / 2
    h = (b - a) / 4
    fml = f(a + h)
    fmr = f(b - h)
    i1 = h / 1.5 * (fa + 4 * fm + fb)
    i2 = h / 3 * (fa + 4 * (fml + fmr) + 2 * fm + fb)
    #i1 = (16 * i2 - i1) / 15
    if ((is + (i1 - i2)) == is) || (m <= a) || (b <= m)
        if ((m <= a) || (b <= m))
            #Interval contains no more machine numbers, required tolerance may not be met
        end
        return i2
    end

    if bottom <= 0
        # println("Number of maximum recursions reached")
        return i2
    end

    return adaptiveSimpsonAux(f, a, m, fa, fml, fm, is, bottom - 1) +
           adaptiveSimpsonAux(f, m, b, fm, fmr, fb, is, bottom - 1)
end
