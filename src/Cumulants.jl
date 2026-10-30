#Cumulants c1, c2 (and c4) of the log-forward from the Taylor expansion of the log characteristic function at 0.
#Models with closed-form cumulants (Heston, Schobel-Zhu, CGMY, ...) specialize computeCumulants.
using TaylorSeries
function computeCumulants(cf::CharFunc{MT,CR}, τ::T; n=4) where {MT,CR,T}
    t = Taylor1(T, n + 1)
    cft = evaluateLogCharFunc(cf, t, τ)
    if !(cft isa Taylor1) || length(cft.coeffs) < n + 1
        throw(ArgumentError(string("cumulants are not available for ", nameof(typeof(cf)), ", use an explicit truncation range [a, b] instead")))
    end
    #println(cft)
    if n == 2
        return (imag(cft.coeffs[2]), -real(cft.coeffs[3]) * 2)
    elseif n == 4
        return (imag(cft.coeffs[2]), -real(cft.coeffs[3]) * 2, real(cft.coeffs[5] * 2 * 3 * 4))
    else
        throw(DomainError(n, "n must be 2 or 4 for cumulants, got $n"))
    end
end
