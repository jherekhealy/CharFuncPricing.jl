export CVCharFunc, CVKind, InitialControlVariance, FullControlVariance, JoshiYangControlVariance

abstract type AbstractTransformation end
struct ALTransformation <: AbstractTransformation
end
transform(::ALTransformation, z::T) where {T} = (1 + z) / (1 - z)
inverseTransform(::ALTransformation, x::T) where {T} = (x - 1) / (x + 1)
transformDerivative(::ALTransformation, z::T) where {T} = 2 / (1 - z)^2
interval(::ALTransformation, ::Type{T}) where {T} = (-one(T), one(T))
isinfTransform(::ALTransformation, z::T) where {T} = z == one(T)
struct LogTransformation{T} <: AbstractTransformation
	c∞::T
end
transform(t::LogTransformation, z::T) where {T} = -log(-z) / t.c∞
inverseTransform(t::LogTransformation, x::T) where {T} = -exp(-x * t.c∞)
transformDerivative(t::LogTransformation, z::T) where {T} = -one(T) / (z * t.c∞)
interval(::LogTransformation{T}, ::Type{T}) where {T} = (-one(T), zero(T))
isinfTransform(::LogTransformation, z::T) where {T} = z == zero(T)
struct IdentityTransformation{T} <: AbstractTransformation
	a::T
	b::T
end
transform(t::IdentityTransformation, z::T) where {T} = z
inverseTransform(t::IdentityTransformation, x::T) where {T} = x
transformDerivative(t::IdentityTransformation, z::T) where {T} = one(T)
interval(t::IdentityTransformation{T}, ::Type{T}) where {T} = (t.a, t.b)
isinfTransform(::IdentityTransformation, z::T) where {T} = false

struct LinearTransformation{T} <: AbstractTransformation
	slope::T
	ordinate::T
end

transform(t::LinearTransformation, z::T) where {T} = t.slope * z + t.ordinate
inverseTransform(t::LinearTransformation, x::T) where {T} = (x - t.ordinate) / t.slope
transformDerivative(t::LinearTransformation, z::T) where {T} = t.slope
interval(t::LinearTransformation{T}, ::Type{T}) where {T} = (t.a, t.b)
isinfTransform(::LinearTransformation, z::T) where {T} = false



struct BlackParams{T}
	σ::T
end

abstract type CVKind end
struct InitialControlVariance <: CVKind end
struct FullControlVariance <: CVKind end
struct JoshiYangControlVariance <: CVKind end
struct AttariControlVariance <: CVKind end

struct CVCharFunc{MAINT, CONTROLT, CR} <: CharFunc{MAINT, CR}
	main::CharFunc{MAINT, CR}
	control::CharFunc{CONTROLT, CR}
end

model(cf::CVCharFunc) = model(cf.main)
oneim(cf::CVCharFunc) = oneim(cf.main)


getControlVariance(::CharFunc{MT}) where {MT} = 0.0

@inline function getControlVariance(
	cf::CVCharFunc{MT, BlackParams{T}})::T where {MT, T}
	return model(cf.control).σ^2
end


function evaluateCharFunc(
	p::CVCharFunc{MAINT, CONTROLT, CR},
	z::CT,
	τ::T) where {T, CR, CT, MAINT, CONTROLT}
	phi = evaluateCharFunc(p.main, z, τ)
	phiB = evaluateCharFunc(p.control, z, τ)
	return phi - phiB
end

function evaluateCharFuncAndDerivative(
	p::CVCharFunc{MAINT, CONTROLT, CR},
	z::CT,
	τ::T,
)::Tuple{CR, CR} where {T, CR, CT, MAINT, CONTROLT}
	phi, phi_d = evaluateCharFuncAndDerivative(p.main, z, τ)
	phiB, phiB_d = evaluateCharFuncAndDerivative(p.control, z, τ)
	return (phi - phiB, phi_d - phiB_d)
end



DefaultCharFunc(params::BlackParams{T}) where {T} = DefaultCharFunc{BlackParams{T}, Complex{T}}(params)

function evaluateLogCharFunc(
	p::CVCharFunc{MAINT, CONTROLT, CR},
	z::CT,
	τ::T) where {T, CR, CT, MAINT, CONTROLT}
	return log(evaluateCharFunc(p, z, τ))
end

function evaluateLogCharFunc(
	cf::CharFunc{BlackParams{T}, CR},
	z::CT,
	τ::T,
) where {T, CR, CT}
	cc1 = oneim(cf)
	p = cf.model
	return -p.σ^2 * τ / 2 * z * (z + cc1)
end

function evaluateLogCharFuncAndDerivative(
	cf::CharFunc{BlackParams{T}, CR},
	z::CT,
	τ::T,
)::Tuple{CR, CR} where {T, CR, CT}
	cc1 = oneim(cf)
	p = cf.model
	phi = -p.σ^2 * τ / 2 * z * (z + cc1)
	phi_d = -p.σ^2 * τ / 2 * (2 * z + cc1)
	return phi, phi_d
end



function computeControlVariance(
	cf::CharFunc,
	τ::T, ::JoshiYangControlVariance,
)::T where {T}
	phid0 = imag(derivative(u -> evaluateCharFunc(cf, u - 1im, τ), eps(T))) #ForwardDiff at zero(T)wrong due to ifelse in charfunc eval
	#println("phi'(-i)=", phid0)
	return 2phid0 / τ
end

function computeControlVariance(
	cf::CharFunc,
	τ::T, ::AttariControlVariance,
)::T where {T}
	phid0 = -imag(derivative(u -> evaluateCharFunc(cf, u, τ), eps(T))) #ForwardDiff at zero(T)wrong due to ifelse in charfunc eval
	return 2phid0 / τ
end

function blackScholes(isCall::Bool, strike, forward, variance, discountDf)
	sign = if isCall
		1
	else
		-1
	end
	sqrtVar = sqrt(variance)
	d1 = log(forward / strike) / sqrtVar + sqrtVar / 2
	d2 = d1 - sqrtVar
	nd1 = normcdf(sign * d1)
	nd2 = normcdf(sign * d2)
	price = discountDf * sign * (forward * nd1 - strike * nd2)
	return price
end
