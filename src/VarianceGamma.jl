export VarianceGammaParams, DefaultCharFunc
struct VarianceGammaParams{T}
    θ::T
    σ::T
    ν::T
end

DefaultCharFunc(params::VarianceGammaParams{Float64}) = DefaultCharFunc{VarianceGammaParams{Float64},Complex{Float64}}(params)

function evaluateCharFunc(cf::CharFunc{VarianceGammaParams{T},CR}, z::CT, τ::T) where {T,CR,CT}
    p = model(cf)
    cc1 = oneim(cf)
    iz = cc1 * z
    return (1 - iz * p.ν * (p.θ + p.σ^2 * iz / 2))^(-τ / p.ν) *exp(iz*τ/p.ν*log(1+0*cc1-p.θ*p.ν-p.σ^2*p.ν/2))
    #return (1/(1 - iz*p.θ*p.ν + p.σ^2 * z^2 * p.ν / 2))^(τ/p.ν)
end

function evaluateLogCharFunc(cf::CharFunc{VarianceGammaParams{T},CR}, z::CT, τ::T) where {T,CR,CT}
    p = model(cf)
    cc1 = oneim(cf)
    iz = cc1 * z
    logValue = (-τ / p.ν) * log(1 - iz * p.ν * (p.θ + p.σ^2 * iz / 2)) +iz*τ/p.ν*log(1+0*cc1-p.θ*p.ν-p.σ^2*p.ν/2)
    #return log(return (1/(1 - iz*p.θ*p.ν + p.σ^2 * z^2 * p.ν / 2))^(τ/p.ν))
    return logValue
end

function computeTruncation(cf::CharFunc{VarianceGammaParams{T1}}, τ::T, tol::T) where {T,T1}
    l = -log(tol) / 2
    c1, c2, c4 = computeCumulants(cf, τ)
    c2 = c2 + sqrt(abs(c4))
    b = c1 + l * sqrt(abs(c2))
    println("truncation b=", b)
    return b
end


asymptoticLogOrder(model::VarianceGammaParams{T}) where {T} = one(T) #not true but ok for SWIFT.