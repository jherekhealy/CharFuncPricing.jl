import Nemo
import Nemo: acb, arb, AcbField,  onei
export NemoCharFunc
#support for Nemo

struct NemoCharFunc{MT} <: CharFunc{MT,acb} #model type, return type (e.g. Complex or acb)
    model::MT
    field::AcbField
end

@inline CharFuncPricing.model(cf::NemoCharFunc) = cf.model
@inline CharFuncPricing.oneim(cf::NemoCharFunc) = onei(cf.field)
@inline CharFuncPricing.const_pi(cf::NemoCharFunc) = real(Nemo.const_pi(cf.field))
@inline Base.expm1(x::acb) = exp(x) - 1 #for Nemo - now in latest Nemo 0.19.2
@inline Base.sign(x::arb) = x < 0 ? -1 : 1
