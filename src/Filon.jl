export priceEuropean, AdaptiveFilonCharFuncPricer, FilonCharFuncPricer



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

#Filon with truncation: very fast, but truncation rule not always robust
struct FilonCharFuncPricer{T}
    τ::T
    b::T
    kcos::Array{T,2}
    vControl::T
    pi::T
    function FilonCharFuncPricer(
        p::CharFunc{MAINT,CR},
        τ::T;
        tTol::T = T(1e-8),
        qTol::T = T(1e-8),
        b::T = Base.zero(T),
        maxRecursionDepth::Int=16
    ) where {MAINT,CR, T}
        if b == 0
            b = computeTruncation(p, τ, T(1e-4))
        end
        mcos = Dict{T,CR}()
        iPure = oneim(p)
        @inline function integrand(x::T)::CR where {T}
            u = x - 0.5 * iPure
            phi = evaluateCharFunc(p, u, τ) 
            denom = 0.25 + x^2
            phi /= denom
            return phi
        end
        @inline function fcos(x::T)::T where {T}
            v = integrand(x)
            mcos[x] = v
            return real(v)
        end
        @inline function fsin(x::T)::T where {T}
            # v = get!(mcos,x , integrand(x))
            local v
            if haskey(mcos, x)
                v = mcos[x]
            else
                v = integrand(x)
                mcos[x] = v
            end
            return imag(v)
        end
        #TODO manual stack management.
        local a = T(0)
        local ic = integrateSimpsonGG(fcos, a,b, qTol, maxRecursionDepth=maxRecursionDepth)
        local is = integrateSimpsonGG(fsin, a,b, qTol, maxRecursionDepth=maxRecursionDepth)
    #    modsim(fcos, a,b, qTol, maxRecursionDepth=16)
    #    modsim(fsin, a,b, qTol, maxRecursionDepth=16)      
   
        #test if truncation is ok
        local sortedKeys
        for tailIteration = 1:24
            sortedKeys = sort!(collect(keys(mcos)))
            n = length(sortedKeys)
            an = sortedKeys[n-2]
            cn = sortedKeys[n-1]
            bn = sortedKeys[n]
            fan = mcos[an]
            fbn = mcos[bn]
            fcn = mcos[cn]
            tailEstimate = simpsonAux(an, bn, real(fan), real(fbn), real(fcn))
            tailEstimateS = simpsonAux(an, bn, imag(fan), imag(fbn), imag(fcn))
            #  println(tailIteration," tail ",tailEstimate," ", tailEstimate/(bn-an)," ", ic, " ",ic*qTol," b=",b," ",bn-an)
            if (abs(tailEstimate) / min(T(1),bn - an) > tTol *  abs(ic)  || abs(tailEstimateS) / min(T(1),bn - an) > tTol *  abs(is)) && tailIteration < 24
                a = b
                b *= 2
                depth = maxRecursionDepth #ceil(Int,maxRecursionDepth/2)
                ic += integrateSimpsonGG(fcos, a, b, qTol, maxRecursionDepth=depth, integralEstimate=ic)
                is += integrateSimpsonGG(fsin, a, b, qTol, maxRecursionDepth=depth, integralEstimate=is)
            else
                break
            end
        end
        kcos = zeros(T, (3, length(mcos)))
        for (i, u) in enumerate(sortedKeys)
            v = mcos[u]
            kcos[:, i] = [u, real(v),imag(v)]
        end

        # for (i, u) in enumerate(sortedKeys)
        #     v = mcos[u]
        #      kcos[:, i*2-1] = [x, real(v),imag(v)]
        # end
        # #add mid points
        # for i = 1:length(mcos)-1
        #     a = @view kcos[:, i*2-1]
        #     b = @view kcos[:, (i+1)*2-1]
        #     x = (a[1] + b[1]) / 2
        #     v = integrand(u)
        #     kcos[:, i*2] = [x, real(v), imag(v)]
        # end

        return new{T}(τ, b, kcos, getControlVariance(p), const_pi(p))
    end
end

function priceEuropean(
    p::FilonCharFuncPricer{T},
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
    local fp::FilonParams
    @inbounds for i = 1:2:nx-2
        if abs(diff - (kcos[1, i+2] - kcos[1, i])) > 1e-8
            fp = FilonParams(freq, kcos[1, i], kcos[1, i+2])
            diff = kcos[1, i+2] - kcos[1, i]
        end
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
