module CharFuncPricing
export DefaultCharFunc, CharFunc

abstract type CharFunc{MT,CR} end

struct DefaultCharFunc{MT,CR} <: CharFunc{MT,CR} #model type, return type (e.g. Complex or acb)
    model::MT
end

@inline model(cf::DefaultCharFunc) = cf.model
@inline oneim(cf::CharFunc) = 1im
@inline oneim(cf::CharFunc{MT,Complex{T}}) where {MT,T} = 1im
@inline oneim(cf::CharFunc{MT,Complex}) where {MT} = 1im
@inline Base.zero(cf::CharFunc{MT,Complex}) where {MT} = Base.zero(Float64)
@inline Base.zero(cf::CharFunc{MT,Complex{BigFloat}}) where {MT} = Base.zero(BigFloat)
@inline const_pi(cf::CharFunc{MT,Complex{T}}) where {MT,T} = T(pi)
@inline const_pi(cf::CharFunc{MT,Complex}) where {MT} = Float64(pi)



include("lambertw.jl")
include("CVCharFunc.jl")
include("Heston.jl")
include("DoubleHeston.jl")
include("SchobelZhu.jl")
include("SVCJ.jl")
include("CGMY.jl")
include("Cos.jl")
include("CosLipton.jl")
include("CosAlpha.jl")
include("quadratures.jl")
include("Flinn.jl")
include("Filon.jl")
include("NemoCharFunc.jl")
include("GaussLobatto.jl")
include("AndersenLake.jl")
include("JoshiYang.jl")
include("Swift.jl")
include("ChebyshevCharFunc.jl")
end
