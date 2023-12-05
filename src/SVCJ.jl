export SVCJParams
#Duffie stochastic volatility with correlated jumps.
struct SVCJParams{T}
    heston::HestonParams{T}
    jumpIntensity::T
    μBar::T
    sjumpVariance::T
    vjumpMean::T
    jumpCorrelation::T
end

DefaultCharFunc(params::SVCJParams{Float64}) = DefaultCharFunc{SVCJParams{Float64},Complex}(params)

function evaluateLogCharFunc(p::CharFunc{SVCJParams{TT},CR}, z::CT, T::TT) where {TT,CR,CT}
    modelj = CharFuncPricing.model(p)

    V0 = modelj.heston.v0
    r = zero(TT) 
    q = zero(TT) 
    kappav = modelj.heston.κ
    rho = modelj.heston.ρ
    thetav = modelj.heston.θ
    sigmav = modelj.heston.σ
    lambda = modelj.jumpIntensity
    mubar = modelj.μBar
    mus = log((1 + modelj.μBar) * (1 - modelj.jumpCorrelation * modelj.vjumpMean)) - modelj.sjumpVariance / 2
    rhoj = modelj.jumpCorrelation
    sigmas = sqrt(modelj.sjumpVariance)
    muv = modelj.vjumpMean
    nubar = thetav
    sigmacy = sigmas
    mucy = mus
    mucv = muv
    rhobar = rho
    zetabar = q
    lambdac = lambda
    # not needed for SVJJ, ie SVCJ
    sigmay = 0.0
    lambday = 0.0
    lambdav = 0.0
    muv = 0.0
    muy = 0.0

    lambdabar = lambday + lambdav + lambdac

    #  thetayfunc = function (c)
    #     exp(muy*c+1/2*sigmay^2*c^2)
    # end

    # thetavfunc = function(c)
    #     1/(1-muv*c)
    # end

    #  thetacfunc = function(c1,c2)
    #      exp(mucy*c1+1/2*sigmacy^2*c1^2) /             (1-mucv*c2-rhoj*mucv*c1)
    # end

    #  thetafunc = function(c1, c2)
    #     (lambday*thetayfunc(c1)+lambday*thetavfunc(c2)+lambdac*thetacfunc(c1,c2))/lambdabar
    # end

    #  mubar = thetafunc(1,0) - 1
    # println(mubar," ",modelj.μBar)

    fy = function (u, tau)
        tau * exp(muy * u +  sigmay^2 / 2 * u .^ 2 )
    end
    fv = function (u, tau)
        a = u .* (1 - u)
        b = sigmav * rhobar * u - kappav
        gamma = sqrt(b .^ 2 + a * sigmav^2)
        if (a == 0)
            return (gamma - b) ./ ((gamma - b) + muv * a) * tau  - (2 * muv ) ./ (sigmav^2 + 2*b * muv ) .* log(1 - ((gamma + b) - muv * a) ./ (2 * gamma) .* (1 - exp(-gamma * tau)))
        end
        return (gamma - b) ./ ((gamma - b) + muv * a) * tau - (2 * muv * a) ./ (gamma .^ 2 - (b - muv * a) .^ 2) .* log(1 - ((gamma + b) - muv * a) ./ (2 * gamma) .* (1 - exp(-gamma * tau)))
    end


    fc = function (u, tau)
        a = u .* (1 - u)
        b = sigmav * rhobar * u - kappav
        c_ = 1 - rhoj * mucv * u
        gamma = sqrt(b .^ 2 + a * sigmav^2)
        if a == 0
            return (gamma - b) ./ ((gamma - b) .* c_ + mucv * a) * tau - (2 * mucv ) ./ ((sigmav .* c_) .^ 2 +2*(b .* c_ * mucv )) .* log(1 - ((gamma + b) .* c_ - mucv * a) / (2 * gamma .* c_) .* (1 - exp(-gamma * tau)))
        end
        d_ = (gamma - b) ./ ((gamma - b) .* c_ + mucv * a) * tau - (2 * mucv * a) ./ ((gamma .* c_) .^ 2 - (b .* c_ - mucv * a) .^ 2) .* log(1 - ((gamma + b) .* c_ - mucv * a) / (2 * gamma .* c_) .* (1 - exp(-gamma * tau)))
        return exp(mucy * u + sigmacy^2 * u .^ 2 / 2) .* d_
    end
    alpha0 = function (tau, u)
        a = u .* (1 - u)
        b = sigmav * rhobar * u - kappav
        gamma = sqrt(b .^ 2 + a * sigmav^2)
        return -r * tau + (r - zetabar) * u * tau - kappav * nubar * (((gamma + b) / sigmav^2) * tau + (2 / sigmav^2) * log(1 - ((gamma + b) ./ (2 * gamma)) .* (1 - exp(-gamma * tau))))
    end


    alphabar = function (tau, u)
        thetainter =  (lambday * fy(u, tau) + lambdav * fv(u, tau) + lambdac * fc(u, tau))/lambdabar
        alpha0(tau, u) - lambdabar * tau * (1 + mubar * u) + lambdabar * thetainter
    end


    betabar = function (tau, u)
        a = u .* (1 - u)
        b = sigmav * rhobar * u - kappav
        gamma = sqrt(b .^ 2 + a * sigmav^2)
        -a .* (1 - exp(-gamma * tau)) ./ (2 * gamma - (gamma + b) .* (1 - exp(-gamma * tau)))
    end


    psi = function (u, v, t, T)
        (alphabar(T - t, u)  + betabar(T - t, u) * v)
    end

    return psi(1im*z,V0,0,T)
end

import DoubleExponentialFormulas: quadde

#Direct implementation of Duffie pricing formulae FIXME remove duplicate code
function priceVanillaDuffie(cf::CharFunc{SVCJParams{TT},CR}, S0, K, T; r=0.0) where {TT,CR}
    y0 = log(S0)
    psichi = function (u, y0, t, T)
       exp(u*y0+evaluateLogCharFunc(cf, -1im*u,T-t))
    end

    Gab = function (a, b, y, T)
        integratefunc = ν -> imag(psichi(complex(a, ν * b), y0, zero(T), T) .* exp(complex(0, -ν * y))) ./ ν
        iv, e = quadde(integratefunc, 0, Inf)
        real(psichi(a, y0, 0, T)/2 - iv / pi)
    end

    aT = -r * T
    bT = 0
    d = 1
   
    GbTplusdminusd = Gab(bT + d, -d, -log(K), T)
    GbTminusd = Gab(bT, -d, -log(K),  T)
    exp(aT) * (GbTplusdminusd - K * GbTminusd)
end

