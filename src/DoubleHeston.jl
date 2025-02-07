export DoubleHestonParams, DoubleHestonCVCharFunc

struct DoubleHestonParams{T}
    heston1::HestonParams{T}
    heston2::HestonParams{T}
end

Base.length(DoubleHestonParams) = 10
Base.iterate(p::DoubleHestonParams, state=1) = 
if state <= 5
    return Base.iterate(p.heston1, state)
elseif state <= 10
    value, newState = Base.iterate(p.heston2, state-5)
    return (value, newState+5)
else 
    nothing
end


DefaultCharFunc(params::DoubleHestonParams{Float64}) =
    DefaultCharFunc{DoubleHestonParams{Float64},Complex{Float64}}(params)

DefaultCharFunc(params::DoubleHestonParams{T}) where {T} =
    DefaultCharFunc{DoubleHestonParams{T},Complex{T}}(params)

cinf(params::DoubleHestonParams{T}, τ::T) where {T} = cinf(params.heston1, τ) + cinf(params.heston2, τ)


@inline evaluateLogCharFunc(
    p::CharFunc{DoubleHestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CT} = evaluateLogCharFuncAL(p, z, τ)


function evaluateLogCharFuncAL(
    cf::CharFunc{DoubleHestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CT}
    p = model(cf)
    cf1 = DefaultCharFunc(p.heston1)
    cf2 = DefaultCharFunc(p.heston2)
    return evaluateLogCharFuncAL(cf1, z, τ) + evaluateLogCharFuncAL(cf2, z, τ)
end

function evaluateLogCharFuncAndDerivative(
    cf::CharFunc{DoubleHestonParams{T},CR},
    z::CT,
    τ::T)::Tuple{CR,CR} where {T,CR,CT} 
    p = model(cf)
    cf1 = DefaultCharFunc(p.heston1)
    cf2 = DefaultCharFunc(p.heston2)
    return evaluateLogCharFuncAndDerivative(cf1, z, τ) .+ evaluateLogCharFuncAndDerivative(cf2, z, τ)
end

function computeCumulants(p::DoubleHestonParams{T}, τ::T) where {T}
    p = model(cf)
    return computeCumulants(p.heston1, τ) + computeCumulants(p.heston2, τ)
end


function DoubleHestonCVCharFunc(heston::CharFunc{DoubleHestonParams{T},CR}, τ, kind::CVKind=InitialControlVariance()) where {T,CR} 
    CVCharFunc{DoubleHestonParams{T},BlackParams{T},CR}(
        heston,
        DefaultCharFunc{BlackParams{T},CR}(
            BlackParams{T}(computeControlVariance(heston, τ, kind))
        )
    )
end

function makeCVCharFunc(heston::CharFunc{DoubleHestonParams{T},CR}, τ, kind::CVKind) where {T,CR}  
    DoubleHestonCVCharFunc(heston, τ, kind)
end


@inline function computeControlVariance(
    cf::CharFunc{DoubleHestonParams{TT}},
    τ::T, kind::Union{FullControlVariance,InitialControlVariance}
)::T where {T,TT}
    p = model(cf)
    cf1 = DefaultCharFunc(p.heston1)
    cf2 = DefaultCharFunc(p.heston2)
    return computeControlVariance(cf1, τ, kind ) + computeControlVariance(cf2, τ, kind ) 
end

#TODO 4th derivative