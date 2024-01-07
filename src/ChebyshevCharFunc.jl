
#struct is subject to changes
struct ChebyshevCharFunc{CF,MT,CR} <:    CharFuncPricing.CharFunc{MT,CR} 
    delegate::CF
    transformation::LinearTransformation{ComplexF64}    
    upperBound::Float64
    coeff::Vector{CR}
end


function ChebyshevCharFunc(delegate::CharFunc{MT,CR} , upperBound::Float64, n::Int, τ::Float64; transformation::LinearTransformation= LinearTransformation(-1.0+0.0im,0.0im)) where {MT,CR} 
    coeff = zeros(CR,n)
    fvalues = zeros(CR,n)
    f = function(x)
        CharFuncPricing.evaluateLogCharFunc(delegate, transform(transformation, x), τ)
    end
    a = 0.0
    b = upperBound
    nodes = (cheb2nodes(Float64, n) .* (b-a) .+ (a+b)) ./ 2
    @. fvalues =f(nodes) #works only with Lewis formula for now.
    cheb2coeff!(coeff, fvalues)
    return ChebyshevCharFunc{typeof(delegate),MT,CR}(delegate, transformation, upperBound, coeff)
end
function CharFuncPricing.evaluateLogCharFunc(p::ChebyshevCharFunc{CF,MT,CR}, z::CT, τ::T) where {T,CT,CF,MT,CR}
    z = inverseTransform(p.transformation,z)    
    if (real(z) > p.upperBound || real(z) < 0.0) 
        println("!!! ",z)
        return one(CR)*(-300)
    end
    cheb2interp(p.coeff, real(z), 0.0, p.upperBound)
end
function CharFuncPricing.getControlVariance(p::ChebyshevCharFunc, τ::T) where {T}
    return getControlVariance(p.delegate)
end

@inline CharFuncPricing.model(cf:: ChebyshevCharFunc{ T} ) where {T} = CharFuncPricing.model(cf.delegate)

##utility methods below

cheb2nodes(T, n::Int) = @. (cos(((1:n) - 1) * T(pi) / (n-1)))

function cheb2coeff!(coeff::AbstractArray{T}, fValues::AbstractArray{TV}) where {T,TV}
    nC = Base.length(fValues)-1
    for sk = 0:nC
        @inbounds sumi = fValues[1] / 2
        @inbounds @simd for i = 2:nC
            sumi += fValues[i] * cos((i - 1) * sk * pi / nC)
        end
        @inbounds sumi += fValues[nC+1] / 2 * cos(sk * pi)
        coeff[sk+1] = 2 * sumi / nC
    end
end

@inline function cheb2interp(coeff::AbstractArray{T}, x::TZ,a::TZ,b::TZ) where {T,TZ}
    y = (2x - a - b)/(b-a)
       return cheb2interp(coeff,y)
   end
   @inline function cheb2interp(coeff::AbstractArray{T}, zck::TZ) where {T,TZ}
       b2 = 0.0
       nC = Base.length(coeff) - 1
       b1 = coeff[nC+1] / 2
       @inbounds @fastmath for sk22 = nC:-1:2
           bd = coeff[sk22] - b2
           b2 = b1
           b1 = 2 * zck * b1 + bd
       end
       b0 = coeff[1] + 2 * zck * b1 - b2
       qck = (b0 - b2) / 2
       qck
   end
