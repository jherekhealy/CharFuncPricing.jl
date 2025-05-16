import StatsBase: mean
export SchobelZhuTSParams

struct SchobelZhuTSParams{T,TTime}
    v0::T
    κ::Vector{T}
    θ::Vector{T}
    ρ::Vector{T}
    σ::Vector{T}
    startTime::Vector{TTime}
end

Base.length(p::SchobelZhuTSParams) = 1 + length(p.κ) * 5


DefaultCharFunc(params::SchobelZhuTSParams{Float64,Float64}) =
    DefaultCharFunc{SchobelZhuTSParams{Float64,Float64},Complex{Float64}}(params)

DefaultCharFunc(params::SchobelZhuTSParams{T,TT}) where {T,TT} =
    DefaultCharFunc{SchobelZhuTSParams{T,TT},Complex{T}}(params)

function SchobelZhuCVCharFunc(sz::CharFunc{SchobelZhuTSParams{T,TT},CR}, τ, kind::CVKind=InitialControlVariance()) where {T,TT,CR}
    CVCharFunc{SchobelZhuTSParams{T,TT},BlackParams{T},CR}(
        heston,
        DefaultCharFunc{BlackParams{T},CR}(
            BlackParams{T}(computeControlVariance(sz, τ, kind)),
        ),
    )
end

function makeCVCharFunc(cz::CharFunc{SchobelZhuTSParams{T,TT},CR}, τ, kind::CVKind) where {T,TT,CR}
    SchobelZhuCVCharFunc(cz, τ, kind)
end


@inline function computeControlVariance(
    cf::CharFunc{SchobelZhuTSParams{T1,T2}},
    τ, kind::InitialControlVariance) where {T1,T2}
    p = model(cf)
    p.v0^2
end



function average(p::SchobelZhuTSParams{T1,T2}, τ::T) where {T1,T2,T}
    lastTimeIndex = searchsortedlast(p.startTime, τ)
    if (p.startTime[lastTimeIndex] - τ >= zero(T) && p.startTime[lastTimeIndex] - τ < sqrt(eps(T)) && lastTimeIndex > 1)
        lastTimeIndex -= 1
    end
    κ = mean(@view(p.κ[1:lastTimeIndex]))
    θ = mean(@view(p.θ[1:lastTimeIndex]))
    ρ = mean(@view(p.ρ[1:lastTimeIndex]))
    σ = mean(@view(p.σ[1:lastTimeIndex]))
    return SchobelZhuParams(p.v0, κ, θ, ρ, σ)
end

cinf(params::SchobelZhuTSParams{T}, τ::T) where {T} = cinf(average(params, τ), τ)


function evaluateLogCharFunc(
    cf::CharFunc{SchobelZhuTSParams{T,TT},CR},
    z::CT,
    τ::T,
) where {T,TT,CR,CT}
    #recursive procedure based on C and D.
    H3 = zero(z)
    H4 = zero(z)
    H5 = zero(z)
    p = model(cf)
    n = length(p.startTime)
    t1 = τ
    for i ∈ n:-1:1
        t0 = p.startTime[i]
        if (t0 < τ)
            dt = min(t1, τ) - t0
            cParams = SchobelZhuParams(p.v0, p.κ[i], p.θ[i], p.ρ[i], p.σ[i])
            H3, H4, H5 = evaluateH3H4H5(cParams, z, H3, H4, H5, dt)
        end
        t1 = t0
    end
    #println(z)
    return p.v0^2 / 2 * H3 + p.v0 * H4 + H5
end

# See "Option Pricing with Piecewise-Constant Parameters, Discrete Jumps and Regime-Switching" Jianwei Zhu, February 3, 2010.
function evaluateH3H4H5(
    p::SchobelZhuParams{T},
    z::CT,
    H30::CR,
    H40::CR,
    H50::CR,
    τ::T,
) where {T,CT,CR}
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = 1im
    iu = cc1 * z
    s12 = -iu / 2 * (iu * (1 - ρ^2) - 1 + 2 * ρ * κ / σ)
    s22 = ρ * κ * θ / σ * iu
    s32 = ρ / (2 * σ) * iu + H30 / 2
    s42 = H40
    γ1 = sqrt(2 * σ^2 * s12 + κ^2)
    γ3 = κ^2 * θ - s22 * σ^2
    γ2 = (κ - 2 * σ^2 * s32) / γ1
    ch1 = cosh(γ1 * τ)
    sh1 = sinh(γ1 * τ)

    γ4 = ch1 + γ2 * sh1
    H3 = (κ - γ1 * (sh1 + γ2 * ch1) / γ4) / σ^2 - 2s32
    H4 = ((κ * θ * γ1 - γ2 * γ3) * (1 - ch1) - (κ * θ * γ1 * γ2 - γ3) * sh1) / (γ4 * γ1 * σ^2) + s42 / γ4
    H5 = -log(γ4) / 2 + ((κ * θ * γ1 - γ2 * γ3)^2 - γ3^2 * (1 - γ2^2)) * sh1 / (γ1^3 * γ4 * σ^2 * 2)
    H5 += (κ * θ * γ1 - γ2 * γ3) * γ3 * (γ4 - 1) / (γ1^3 * σ^2 * γ4) + τ * (κ * γ1^2 * (σ^2 - κ * θ^2) + γ3^2) / (γ1^2 * σ^2 * 2)
    H5 += s42 * (γ3 * (γ4 - 1) + (κ * θ * γ1 + σ^2 * γ1 * s42 / 2 + γ2 * γ3) * sh1) / (γ1^2 * γ4) - s32 * σ^2 * τ
    H5 += H50
    #return -s32 * v0^2 - iu * ρ * σ * τ / 2 + H3 / 2 * v0^2 + H4 * v0 + H5
    return H3, H4, H5
end
