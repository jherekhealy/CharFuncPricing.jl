export CGMYParams, evaluateCharFunc, evaluateLogCharFunc, computeCumulants
import SpecialFunctions: gamma
struct CGMYParams{T}
    C::T
    G::T
    M::T
    Y::T
    σ::T
    ω::T    
 CGMYParams(C::T,G::T,M::T,Y::T) where {T} = new{T}(C,G,M,Y,zero(T),driftCorrection(C,G,M,Y))
 CGMYParams(C::T,G::T,M::T,Y::T,σ::T) where {T} = new{T}(C,G,M,Y,σ,driftCorrection(C,G,M,Y))
end

DefaultCharFunc(params::CGMYParams{Float64}) = DefaultCharFunc{CGMYParams{Float64},Complex}(params)

@inline function evaluateCharFunc(p::CharFunc{CGMYParams{T},CR}, z::CT, τ::T)::CR where {T,CR,CT}
    return exp(evaluateLogCharFunc(p, z, τ))
end


function evaluateLogCharFunc(cf::CharFunc{CGMYParams{T},CR}, z::CT, τ::T)::CR where {T,CR,CT}
    p = model(cf)
    σ = p.σ
    cc1 = oneim(cf)
    iz = cc1*z
    q = p.ω
    return iz*q*τ-(z*σ)^2*τ/2 + p.C*τ*gamma(-p.Y)*( (p.M-iz)^p.Y - p.M^p.Y + (p.G+iz)^p.Y - p.G^p.Y )    
end

function computeCumulants(cf::CharFunc{CGMYParams{T},CR}, τ::T) where {T,CR}
    return computeCumulants(model(cf),τ)
end

function computeCumulants(p::CGMYParams{T}, τ::T) where {T}
    c1 = p.C*τ*gamma(1-p.Y)*(p.M^(p.Y-1)-p.G^(p.Y-1)) #+ driftCorrection(p)*τ
    c2 = p.σ^2 * τ + p.C * τ*gamma(2-p.Y)*(p.M^(p.Y-2)+p.G^(p.Y-2))
    c4 = p.C* τ*gamma(4-p.Y)*(p.M^(p.Y-4)+p.G^(p.Y-4))
    println(c1, " ",c2," ", c4)
    return c1, c2, c4
end

driftCorrection(p::CGMYParams{T}) where {T} =  driftCorrection(p.C,p.G,p.M,p.Y)
driftCorrection(C::T,G::T,M::T,Y::T) where {T} = -C* gamma(-Y)*((M-one(T))^Y-M^(Y) + (G+one(T))^Y - G^Y)