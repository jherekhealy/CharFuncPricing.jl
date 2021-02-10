import Nemo: acb, arb, AcbField, const_pi, onei
#support for Nemo

@inline CharFuncPricing.pi(cf::CharFunc{MT, acb, AcbField}) where {MT} = real(const_pi(field(cf)))
@inline CharFuncPricing.oneim(cf::CharFunc{MT, acb, AcbField}) where {MT} = onei(field(cf))
@inline Base.expm1(x::acb) = exp(x) - 1 #for Nemo - now in latest Nemo 0.19.2
Base.sign(x::arb) = x < 0 ? -1 : 1
@inline Base.zero(cf::CharFunc{MT, acb, AcbField}) where {MT} = real(field(cf)(0))
