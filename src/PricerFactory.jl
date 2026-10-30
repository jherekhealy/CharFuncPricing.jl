export makeCharFuncPricer, pricerMethod

#Prices for many strikes. Pricers may specialize this method when they can share work across strikes (e.g. CosFFTCharFuncPricer).
function priceEuropean(p, isCalls::AbstractArray{Bool}, strikes::AbstractArray, forward, τ, discountDf)
    return map((isCall, strike) -> priceEuropean(p, isCall, strike, forward, τ, discountDf), isCalls, strikes)
end

#accepts symbols or strings, with the older hyphenated names, e.g. "Joshi-Yang" for :JoshiYang
pricerMethod(method) = Val(Symbol(replace(string(method), "-" => "")))

makeCharFuncPricer(method::Union{Symbol,AbstractString}, cf, t) = makeCharFuncPricer(pricerMethod(method), cf, t)

function makeCharFuncPricer(::Val{:AndersenLake}, cf, t)
    ALCharFuncPricer(cf, n=64)
end

function makeCharFuncPricer(::Val{:JoshiYang}, cf, t)
    JoshiYangCharFuncPricer(cf, t, n=64)
end
function makeCharFuncPricer(::Val{:Flinn}, cf, t)
    FlinnCharFuncPricer(cf, t, tTol=1e-4, qTol=1e-8)
end

function makeCharFuncPricer(::Val{:FlinnTransformed}, cf, t)
    AdaptiveFlinnCharFuncPricer(cf, t)
end
function makeCharFuncPricer(::Val{:Swift}, cf, t)
    m, _ = findSwiftScaling(cf, t)
    makeSwiftCharFuncPricer(cf, t, m, 3)
end
function makeCharFuncPricer(::Val{:Cos}, cf, t)
    makeCosCharFuncPricer(cf, t, 256, 16)
end
function makeCharFuncPricer(::Val{:Cos128}, cf, t)
    makeCosCharFuncPricer(cf, t, 128, 12)
end
function makeCharFuncPricer(::Val{:CosJunike}, cf, t)
    makeCosCharFuncPricer(cf, t, tol=1e-4, maxM=4096)
end
function makeCharFuncPricer(::Val{:FCC}, cf, t)
    FCCCharFuncPricer(cf, t, qTol=1e-8)
end
