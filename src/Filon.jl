export priceEuropean, AdaptiveFilonCharFuncPricer

abstract type AbstractTransformation end
struct ALTransformation <: AbstractTransformation
end
transform(::ALTransformation, z::T) where {T} = (1 + z) / (1 - z)
inverseTransform(::ALTransformation, x::T) where {T} = (x - 1) / (x + 1)
transformDerivative(::ALTransformation, z::T) where {T} = 2 / (1 - z)^2
interval(::ALTransformation, ::Type{T}) where {T} = (-one(T),one(T))
isinfTransform(::ALTransformation, z::T) where {T} =  z == one(T)
struct LogTransformation{T} <: AbstractTransformation
    c∞ ::T
end
transform(t::LogTransformation, z::T) where {T} = -log(-z)/t.c∞
inverseTransform(t::LogTransformation, x::T) where {T} = -exp(-x * t.c∞)
transformDerivative(t::LogTransformation, z::T) where {T} = -one(T)/(z*t.c∞)
interval(::LogTransformation{T}, ::Type{T}) where {T} = (-one(T),zero(T))
isinfTransform(::LogTransformation, z::T) where {T} =  z == zero(T)
struct IdentityTransformation{T} <: AbstractTransformation
    a ::T
    b ::T
end
transform(t::IdentityTransformation, z::T) where {T} = z
inverseTransform(t::IdentityTransformation, x::T) where {T} = x
transformDerivative(t::IdentityTransformation, z::T) where {T} = one(T)
interval(t::IdentityTransformation{T}, ::Type{T}) where {T} = (t.a, t.b)
isinfTransform(::IdentityTransformation, z::T) where {T} =  false
# struct CotTransformation <: AbstractTransformation
# end
# transform(t::CotTransformation, z::T) where {T} = cot(z/2)^2
# inverseTransform(t::CotTransformation, x::T) where {T} = acot(sqrt(x))*2
# transformDerivative(t::CotTransformation, z::T) where {T} = 2*sin(z)/(1-cos(z))^2
# interval(t::CotTransformation, ::Type{T}) where {T} = (zero(T),T(pi))
# isinfTransform(::CotTransformation, z::T) where {T} =  z == zero(T) 
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
        myTrans = ALTransformation(),
        maxRecursionDepth=16
    ) where {T}
        mcos = Dict{T,Tuple{T,T}}()
        iPure = oneim(p)
        ###### PROBLEMS - if we go to one interval (first one) deep, calc stops for all. But we really have a badly distributed set of points
        #                 the transform is good to find interval truncation, but not good for the points if transform is too different from identity.        
        #myTrans = LogTransformation(cinf(model(p), τ)/2) # without /2 fsin not well behaved near z = 0 (x=infty).
        @inline function integrand(z::T)::Tuple{T,T} where {T}
            if isinfTransform(myTrans, z)
                return zero(T), zero(T)
            end
            x = transform(myTrans,z)
            u = x - iPure / 2
            phi = evaluateCharFunc(p, u, τ)
            denom = (x^2 + 1 / 4)
            phi /= denom
            factor = transformDerivative(myTrans,z)
            phi *= factor
            return real(phi), imag(phi)
        end

        @inline function fcos(x::T)::T where {T}
            v = integrand(x)
            if  !isinfTransform(myTrans, x)
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
                if !isinfTransform(myTrans,x)
                    mcos[x] = v
                end
            end
            return v[2]
        end
        (a,b) = interval(myTrans, T)
        integrateSimpsonGG(fcos, a,b, qTol, maxRecursionDepth=maxRecursionDepth)
        integrateSimpsonGG(fsin, a,b, qTol, maxRecursionDepth=maxRecursionDepth)
    #    modsim(fcos, a,b, qTol, maxRecursionDepth=16)
    #    modsim(fsin, a,b, qTol, maxRecursionDepth=16)
       
        sortedKeys = sort!(collect(keys(mcos)))
        kcos = zeros(T, (3, length(mcos) * 2 - 1))
        #convert to original coordinate
        if (transformDerivative(myTrans, (a+b)/2) < zero(T))
            sortedKeys = reverse(sortedKeys)
        end
        for (i, u) in enumerate(sortedKeys)
            v = mcos[u]
            x = transform(myTrans,u)
            factor = transformDerivative(myTrans, u)
            kcos[:, i*2-1] = [x, v[1] / factor, v[2] / factor]
        end
        #add mid points
        for i = 1:length(mcos)-1
            a = @view kcos[:, i*2-1]
            b = @view kcos[:, (i+1)*2-1]
            x = (a[1] + b[1]) / 2
            u = inverseTransform(myTrans, x) #(1-u)*x = 1+u => u(1+x) = x-1 => u=(x-1)/(x+1)
            v = integrand(u)
            factor = transformDerivative(myTrans, u)
            kcos[:, i*2] = [x, v[1] / factor, v[2] / factor]
        end
        #end interval = kcos[1,end]
        return new{T}(τ, kcos, getControlVariance(p), const_pi(p))
    end

end
Base.broadcastable(p::AdaptiveFilonCharFuncPricer) = Ref(p)


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
