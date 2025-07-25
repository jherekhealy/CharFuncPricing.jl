#Implementation of "FOURIER TRANSFORMS, OPTION PRICING AND CONTROLS" by MARK JOSHI AND CHAO YANG
# bad for long maturities e.g. (30.0, 100.0001, HestonParams{Float64}(1.0, 0.5, 1.0, 0.95, 1.0)) as cv explodes

using FastGaussQuadrature
import ForwardDiff: derivative
export JoshiYangCharFuncPricer

struct JoshiYangCharFuncPricer{T, CR}
	τ::T
	w::Array{T, 1}
	x::Array{T, 1}
	phi::Array{CR, 1}
	vControl::T
	pi::T
	function JoshiYangCharFuncPricer(
		cf::CharFunc{MAINT, CR},
		τ::T;
		n = 64,
	) where {MAINT, CR, T}
		blackVariance = abs(computeControlVariance(cf, τ, JoshiYangControlVariance()))
		phicv = makeCVCharFunc(cf, blackVariance)

		x, w = gausslaguerre(n)
		iPure = oneim(cf)
		@inbounds for i ∈ eachindex(w)
			if w[i] != 0 && x[i] < 600
				w[i] *= exp(x[i])
			else
				w[i] = 0.0
			end ##else overflow
		end
		phi = @. (evaluateCharFunc(phicv, x - iPure, τ) / (x * (x - iPure)))

		return new{T, CR}(τ, w, x, phi, blackVariance, const_pi(cf))
	end

	function JoshiYangCharFuncPricer(
		cf::CharFunc{MAINT, CR},
		τ::T,
		x::AbstractArray{T},
		w::AbstractArray{T},
	) where {MAINT, CR, T}
		blackVariance = abs(computeControlVariance(cf, τ, JoshiYangControlVariance()))
		phicv = makeCVCharFunc(cf, blackVariance)

		iPure = oneim(cf)
		phi = @. (evaluateCharFunc(phicv, x - iPure, τ) / (x * (x - iPure)))

		return new{T, CR}(τ, w, x, phi, blackVariance, const_pi(cf))
	end

end
Base.broadcastable(p::JoshiYangCharFuncPricer) = Ref(p)
Base.broadcastable(p::CVCharFunc) = Ref(p)
Base.broadcastable(p::DefaultCharFunc) = Ref(p)


makeCVCharFunc(cf::CharFunc{MAINT, CR}, blackVariance::T) where {MAINT, CR, T} =
	if blackVariance == zero(T)
		cf
	else
		blackCf = DefaultCharFunc{BlackParams{T}, CR}(
			BlackParams{T}(sqrt(blackVariance)),
		)
		CVCharFunc(cf, blackCf)
	end



function priceEuropean(
	p::JoshiYangCharFuncPricer{T, CR},
	isCall::Bool,
	strike::T,
	forward::T,
	τ::T,
	discountDf::T,
) where {T, CR}
	if τ != p.τ
		throw(DomainError(τ, string("maturity is different from pricer maturity ", p.τ)))
	end
	#use forward or strike ?   freq= -freq and *strike for strike.

	#psih(xi= exp(i xi log(forward))

	freq = log(forward / strike)
	integral = sum([wi * real(phii * exp(xi * freq * 1im)) for (phii, xi, wi) ∈ zip(p.phi, p.x, p.w)])
	callPrice = -discountDf / p.pi * integral * forward
	if p.vControl != 0
		#use control variate
		blackPrice = blackScholes(true, strike, forward, p.vControl * p.τ, discountDf)

		#println("bs ",price," ",callPrice)
		callPrice += blackPrice
	else
		# println("integ ",callPrice)
		callPrice += discountDf * forward / 2
	end
	if isCall
		return callPrice
	end
	putPrice = callPrice - (forward - strike) * discountDf
	return putPrice
end
