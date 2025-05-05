import StatsBase:mean
export HestonTSParams

struct HestonTSParams{T, TTime}
	v0::T
	κ::Vector{T}
	θ::Vector{T}
	ρ::Vector{T}
	σ::Vector{T}
	startTime::Vector{TTime}
end

Base.length(p::HestonTSParams) = 1 + length(p.κ) * 5


DefaultCharFunc(params::HestonTSParams{Float64, Float64}) =
	DefaultCharFunc{HestonTSParams{Float64, Float64}, Complex{Float64}}(params)

DefaultCharFunc(params::HestonTSParams{T, TT}) where {T, TT} =
	DefaultCharFunc{HestonTSParams{T, TT}, Complex{T}}(params)

function HestonCVCharFunc(heston::CharFunc{HestonTSParams{T, TT}, CR}, τ, kind::CVKind = InitialControlVariance()) where {T, TT, CR}
	CVCharFunc{HestonTSParams{T, TT}, BlackParams{T}, CR}(
		heston,
		DefaultCharFunc{BlackParams{T}, CR}(
			BlackParams{T}(computeControlVariance(heston, τ, kind)),
		),
	)
end

function makeCVCharFunc(heston::CharFunc{HestonTSParams{T, TT}, CR}, τ, kind::CVKind) where {T, TT, CR}
	HestonCVCharFunc(heston, τ, kind)
end


@inline function computeControlVariance(
	cf::CharFunc{HestonTSParams{T1, T2}},
	τ, kind::InitialControlVariance) where {T1, T2}
	p = model(cf)
	p.v0
end


@inline function computeControlVariance(
    cf::CharFunc{HestonTSParams{T1,T2}},
    τ::T, kind::FullControlVariance
)::T where {T,T1,T2}
    p = model(cf)
    averageParams = average(p, τ)
    ektm = 1 - exp(-averageParams.κ * τ)
    return ektm * (p.v0 - averageParams.θ) / (averageParams.κ * τ) + averageParams.θ
end

function average(p::HestonTSParams{T1,T2},τ::T) where {T1,T2,T}
    lastTimeIndex = searchsortedlast(p.startTime,τ)
    if (p.startTime[lastTimeIndex] - τ >= zero(T) &&  p.startTime[lastTimeIndex] - τ < sqrt(eps(T)) &&  lastTimeIndex > 1)
        lastTimeIndex -=1
    end
    κ = mean(@view(p.κ[1:lastTimeIndex]))
    θ = mean(@view(p.θ[1:lastTimeIndex]))
    ρ = mean(@view(p.ρ[1:lastTimeIndex]))
    σ = mean(@view(p.σ[1:lastTimeIndex]))
    return HestonParams(p.v0, κ, θ, ρ, σ)
end

cinf(params::HestonTSParams{T}, τ::T) where {T} = cinf(average(params, τ),τ)


function evaluateLogCharFunc(
	cf::CharFunc{HestonTSParams{T, TT}, CR},
	z::CT,
	τ::T,
) where {T, TT, CR, CT}
	#recursive procedure based on C and D.
	D = zero(z)
	C = zero(z)
	p = model(cf)
	n = length(p.startTime)
    t1 = τ
	for i ∈ n:-1:1
		t0 = p.startTime[i]
		if (t0 < τ)
			dt = min(t1,τ) - t0
			cParams = HestonParams(p.v0, p.κ[i], p.θ[i], p.ρ[i], p.σ[i])
			C, D = evaluateCD(cParams, z, C, D, dt)
		end
        t1 = t0
	end
    #println(z)
	return p.v0 * D + C
end

function evaluateCD(
	p::HestonParams{T},
	z::CT,
	C0::CR,
	D0::CR,
	τ::T,
) where {T, CT,CR}
	κ = p.κ
	θ = p.θ
	ρ = p.ρ
	σ = p.σ
	cc1 = 1im
	α = -(z * (z + cc1)) * σ^2
	β = κ - σ * ρ * z * cc1
	d = sqrt(β^2 - α)
	g = (β - d) / (β + d)
	gtilde = (β - d - D0 * σ^2) / (β + d - D0 * σ^2)
	em = exp(-d * τ)
	D = (β + d) / σ^2 * (g - gtilde * em) / (1 - gtilde * em)
	C = κ * θ / σ^2 * (-2 * log((1 - gtilde * em) / (1 - gtilde)) + (β - d) * τ) + C0
	return (C, D)
	# C + D * v0
end
