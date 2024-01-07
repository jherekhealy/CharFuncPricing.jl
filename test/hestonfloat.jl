using CharFuncPricing, Test, Printf
using StatsBase

@testset "HestonLong" begin
τ=10.0
forward= 10000.0
params=HestonParams(1.0, 2.0, 0.0025, 0.5, 0.1)
strikes = [100.0001, 101.0, 110.0, 200.0, 1000.0, 10000.0]
cf = DefaultCharFunc(params)
quadratureGL = ModlobQuadrature(1e-13)
m, stol = CharFuncPricing.findSwiftScaling(cf,τ)
pricers = [ALCharFuncPricer(cf, quadratureGL), makeCosCharFuncPricer(cf, τ, 200, 8),
 makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10),
 FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, n=200),
 JoshiYangCharFuncPricer(cf,τ,n=16),CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, m, 4)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)","Joshi-Yang (n=64)","Swift"]
isCall = false
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
end
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
#@test isapprox(prices[1]-refPrices[1],1.97, atol=1e-2) #bad price for Flinn
@test isapprox(refPrices[1],prices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[6], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[7], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-15)


τ=10.0
forward = 100.0
params = HestonParams(0.0001, 0.1, 0.25, 0.95, 3.0)
strike = 100.0001
cf = DefaultCharFunc(params)
m, stol = CharFuncPricing.findSwiftScaling(cf,τ)
pricers = [ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
makeCosCharFuncPricer(cf, τ, 200, 8), makeCosCharFuncPricer(cf, τ, 1000, 8),
 FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10),ALCharFuncPricer(cf, n=200),
 JoshiYangCharFuncPricer(cf,τ, n=64),
 CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, m, 4)]
 pricerNames = ["Reference",  "Cos-Adaptive","Cos (M=200)","Cos (M=1000)",
  "Flinn-Truncated (1e-10)",  "Andersen-Lake (n=200)","Joshi-Yang","Swift"]

isCall = false
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
end
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1]-refPrices[1],0.27, atol=1e-2) #bad price for Cos
prices = map(x -> priceEuropean(pricers[2], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-4) # bad accuracy for adaptive Cos
prices = map(x -> priceEuropean(pricers[5], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-7)

#10.0 0.0001 0.01 1.0 -0.95 3.0 10000.0 9900.0 9882.675792888303
params = HestonParams(0.0001, 0.01, 1.0, -0.95, 3.0)
cf = DefaultCharFunc(params)
m, stol = CharFuncPricing.findSwiftScaling(cf,τ)
pricers=[ALCharFuncPricer(cf, quadratureGL), makeCosCharFuncPricer(cf, τ, 8),
 FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10), ALCharFuncPricer(cf, n=200), CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, m, 4)]
 pricerNames=["Reference", "Cos-Adaptive 8", "Flinn Trunc 1e-10", "AndersenLake","Swift"]
isCall = false
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
end
prices = map(x -> priceEuropean(pricers[2], isCall, x, forward, τ, 1.0), strikes)
@test isless(1e-4, abs(prices[1]-refPrices[1])) #bad price for Cos
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)

τ=10.0
strike=101.0
params=HestonParams{Float64}(0.0001, 0.01, 0.0001, 0.95, 3.0)
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, quadratureGL), makeCosCharFuncPricer(cf, τ, 200, 8), makeCosCharFuncPricer(cf, τ, 8, tol=1e-8),
FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10), FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10),
 AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[2], isCall, x, forward, τ, 1.0), strikes)
@test isless(6.28, abs(prices[1]-refPrices[1])) #bad price for Cos
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-4)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
#@test isless(1.86, abs(prices[1]-refPrices[1])) #bad price for Flinn-truncated
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[6], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[7], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)

τ=10.0
strike=101.0
params= HestonParams{Float64}(1.0, 0.01, 0.25, 0.95, 3.0)
cf = DefaultCharFunc(params)
pricers= [ ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
FlinnCharFuncPricer(HestonCVCharFunc(cf,τ), τ, tTol = 1e-10, qTol=1e-10), 
ALCharFuncPricer(cf, n=200), JoshiYangCharFuncPricer(cf, τ, n=128)]
pricerNames=["Reference","Cos Adapt","Flinn trunc", "AL", "Joshi"]
isCall = false
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)

#(30.0, 100.0001, HestonParams{Float64}(1.0, 0.1, 1.0, 0.95, 0.1))
τ=30.0
params =  HestonParams{Float64}(1.0, 0.1, 1.0, 0.95, 0.1)
cf = DefaultCharFunc(params)
pricers= [ ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
FlinnCharFuncPricer(HestonCVCharFunc(cf,τ), τ, tTol = 1e-10, qTol=1e-10), ALCharFuncPricer(cf, n=200)]
pricerNames=["Reference","Cos Adapt","Flinn trunc", "AL"]
isCall = false
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)


τ=0.0025
params = HestonParams{Float64}(0.0001, 0.01, 0.0001, -0.95, 3.0)
cf = DefaultCharFunc(params)
pricers= [ ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
FlinnCharFuncPricer(HestonCVCharFunc(cf, τ), τ, tTol = 1e-10, qTol=1e-10), ALCharFuncPricer(cf, n=10000)]
pricerNames=["Reference","Cos Adapt","Flinn trunc", "AL"]
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-10)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-15)


#flinn breaks strongly with or without cv. Also due to bad truncatino range (better if we use short mat criteria as well)
τ=0.5
strike = 100.0
params = HestonParams{Float64}(0.0025, 2.0, 0.25, 0.5, 0.0001)
cf = DefaultCharFunc(params)
pricers= [ ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
FlinnCharFuncPricer(HestonCVCharFunc(cf, τ), τ, tTol = 1e-10, qTol=1e-10), ALCharFuncPricer(cf, n=10000)]
pricerNames=["Reference","Cos Adapt","Flinn trunc", "AL"]
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-13)



#DE adaptive largest error abs
τ= 0.1
strike= 101.0
forward=100.0
params = HestonParams{Float64}(0.0001, 2.0, 0.04, -0.5, 3.0)
cf = DefaultCharFunc(params)
pricers= [ ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
FlinnCharFuncPricer(HestonCVCharFunc(cf, τ), τ, tTol = 1e-10, qTol=1e-10), ALCharFuncPricer(cf, n=10000)]
pricerNames=["Reference","Cos Adapt","Flinn trunc", "AL"]
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-9)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-14)


#Interestingly, AL-10000 is likley wrong as AL-0 with tol=1e-20 = TFLinn price, strike = 101 & 200.
τ=30.0
forward = 100.0
params = HestonParams{Float64}(0.0025, 0.1, 0.0001, 0.1, 3.0)
cf = DefaultCharFunc(params)
pricers= [ ALCharFuncPricer(cf, quadratureGL),  makeCosCharFuncPricer(cf, τ, 8),
FlinnCharFuncPricer(HestonCVCharFunc(cf, τ), τ, tTol = 1e-10, qTol=1e-10), ALCharFuncPricer(cf, n=10000)]
pricerNames=["Reference","Cos Adapt","Flinn trunc", "AL"]
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-12)


#Largest error on AL-1000 vs 2000
τ=10.0
forward = 100.0
params = HestonParams{Float64}(0.0001, 0.01, 0.0001, 0.95, 1.0)
#(10.0, 100.0001, HestonParams{Float64}(0.0001, 0.01, 0.0001, 0.95, 1.0))
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, quadratureGL), makeCosCharFuncPricer(cf, τ, 200, 8),
makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10),
 FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, τ, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, τ, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
prices = map(x -> priceEuropean(pricers[2], isCall, x, forward, τ, 1.0), strikes)
@test isless(0.95, abs(prices[1]-refPrices[1])) #bad price for Cos
prices = map(x -> priceEuropean(pricers[3], isCall, x, forward, τ, 1.0), strikes)
@test isless(2e-4, abs(prices[1]-refPrices[1])) #bad price for Cos
prices = map(x -> priceEuropean(pricers[4], isCall, x, forward, τ, 1.0), strikes)
#@test isless(5e-2, abs(prices[1]-refPrices[1])) #bad price for Flinn
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[5], isCall, x, forward, τ, 1.0), strikes)
#@test isless(5e-2, abs(prices[1]-refPrices[1])) #bad price for Flinn
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[6], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-7)
prices = map(x -> priceEuropean(pricers[7], isCall, x, forward, τ, 1.0), strikes)
@test isapprox(prices[1], refPrices[1], atol=1e-8)

#largest error n=200
τ=2.0
strike = 100.0
forward = 200.0
forwards = [100.0,100.0001, 101.0, 110.0, 200.0, 1000.0, 10000.0]
params = HestonParams{Float64}(0.0001, 2.0, 0.0001, 0.95, 1.0)
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, quadratureGL), makeCosCharFuncPricer(cf, τ, 200, 8),
 makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10),
  FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, strike, x, τ, 1.0),forwards)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, strike, x, τ, 1.0), forwards)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
end
prices = map(x -> priceEuropean(pricers[7], isCall, strike, x, τ, 1.0), forwards)
@test isapprox(prices[5], refPrices[5], atol=1.2e-2)
@test isless(0.01, prices[5]/ refPrices[5])

end


@testset "HestonFlinn" begin
    v0 = 0.0718
    ρ = -0.352
    σ = 0.582
    θ = 0.0762
    κ = 1.542
    params = HestonParams(v0, κ, θ, ρ, σ)

    tol = 1e-8
    τ = 0.1
    forward = 6317.80
    strike = 6317.80
    cf = DefaultCharFunc(params)
    refPricer = makeCosCharFuncPricer(cf, τ, 2048 * 4, 32)

    pricer = FlinnCharFuncPricer(cf, τ, qTol = tol, tTol = tol)
    refPrice = priceEuropean(refPricer, false, strike, forward, τ, 1.0)
    price = priceEuropean(pricer, false, strike, forward, τ, 1.0)
    @test isapprox(refPrice, price, atol = 5e-4)
#    @test isapprox(209.820637, price, atol = 1e-5)

    ccf = HestonCVCharFunc(cf, τ)
    pricer = FlinnCharFuncPricer(ccf, τ, qTol = tol, tTol = tol)
    price = priceEuropean(pricer, false, strike, forward, τ, 1.0)
    @test isapprox(refPrice, price, atol = 1e-6)

    #Make sure BigFloat works
    params = HestonParams(BigFloat("0.0718"), BigFloat("1.542"), BigFloat("0.0762"),BigFloat("-0.352"),BigFloat("0.582"))
    cf = DefaultCharFunc{HestonParams{BigFloat},Complex{BigFloat}}(params)
    τ = BigFloat("0.1")
    strike = BigFloat("6317.80")
    forward = BigFloat("6317.80")
    pricer = AdaptiveFlinnCharFuncPricer(cf, τ, qTol = BigFloat(1e-24))
    priceB = priceEuropean(pricer, false, strike, forward, τ, BigFloat(1.0))
    @test isapprox(refPrice, Float64(priceB), atol = 1e-12)
    pricer = FlinnCharFuncPricer(cf, τ, qTol = BigFloat(1e-24), tTol = BigFloat(1e-24))
    priceC = priceEuropean(pricer, false, strike, forward, τ, BigFloat(1.0))
    @test isapprox(priceB, priceC, atol=BigFloat(1e-18))
end


@testset "HestonLogContinuity" begin
    v0 = 0.08
    ρ = -0.8
    σ = 0.25
    θ = 0.1
    κ = 3.0
    #params from Cui et al for Heston.
    params = HestonParams(v0, κ, θ, ρ, σ)
    #with the original formula from SZ, there is a discontinuity in the log.
    cf = DefaultCharFunc(params)
    sz = CharFuncPricing.evaluateLogCharFuncCui(cf, 2.0 + 0.0im, 15.0)
    lk = CharFuncPricing.evaluateLogCharFunc(cf, 2.0 + 0.0im, 15.0)
    @test isapprox(lk, sz, atol = 1e-13)
    for u = 0.1:0.1:4.0
        lk = CharFuncPricing.evaluateLogCharFunc(cf, u + 0.0im, 15.0)
        sz = CharFuncPricing.evaluateLogCharFuncCui(cf, u + 0.0im, 15.0)
        @test isapprox(lk, sz, atol = 1e-13)
    end
    v0 = 0.010201
    ρ = -0.7
    σ = 0.61
    κ = 6.21
    θ = 0.019
    α = 3.35861
    params = HestonParams{Float64}(sqrt(v0), κ, sqrt(θ), ρ, σ)
    cf = DefaultCharFunc(params)
    sz = CharFuncPricing.evaluateLogCharFunc(cf, 10.0 - (α + 1) * 1im, 10.0)
    lk = CharFuncPricing.evaluateLogCharFuncCui(cf, 10.0 - (α + 1) * 1im, 10.0)
    @test isapprox(lk, sz, atol = 1e-13)
end

@testset "LF2y" begin
    r = 0.0
    q = 0.0
    κ = 1.0
    θ = 0.1
    σ = 1.0
    ρ = -0.9
    v0 = 0.1
    τs = [2.0/52, 2.0/12, 2.0]
    spot = 1.0

    τ = τs[3]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    m = 200
    l = 8
    refPricer = makeCosCharFuncPricer(cf, τ, 2048 * 8, 32)
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    refPrices = Vector{Float64}(undef, 0)
    prices = Vector{Float64}(undef, 0)
    for strike = 0.4:0.012:1.6
        refPrice = priceEuropean(refPricer, false, strike, spot,τ, df)
        price = priceEuropean(pricer, false, strike, spot, τ, df)
        push!(refPrices, refPrice)
        push!(prices, price)
    end
    rmse = rmsd(prices, refPrices)
    println("RMSE ", rmse)
    @test isless(rmse, 1e-4)
    @test isless(1e-5, rmse)

    for τ in τs
        local refPrices = zeros(length(collect(0.4:0.012:1.6)))
        local prices =  zeros(length(collect(0.4:0.012:1.6)))
        measuredTime = @elapsed begin
        pricer = makeCosCharFuncPricer(cf, τ, 2048 * 8, 32)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot,  τ, df)
        refPrices[i] = price
    end
    end
    println("Ref ", measuredTime)

    measuredTime = @elapsed begin
        pricer = makeCosCharFuncPricer(cf, τ, 200,8)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, τ, df)
        prices[i] = price
    end
    end
    rmse = rmsd(prices, refPrices)
    println("Cos ", rmse, " ",measuredTime)

    measuredTime = @elapsed begin
        pricer = FlinnCharFuncPricer(cf, τ)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot,  τ, df)
        prices[i] = price
    end
    end
    rmse = rmsd(prices, refPrices)
    println("Truncated Flinn ", rmse, " ",measuredTime)

    measuredTime = @elapsed begin
        pricer = AdaptiveFlinnCharFuncPricer(cf, τ)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, τ, df)
        prices[i] = price
    end
    end
    rmse = rmsd(prices, refPrices)
    println("Transformed Flinn ", rmse, " ",measuredTime)

quadrature = TanhSinhQuadrature(200, eps())
    measuredTime = @elapsed begin
        pricer = ALCharFuncPricer(cf, quadrature)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, τ, df)
        prices[i] = price
        end
    end
    rmse = rmsd(prices, refPrices)
    println("Andersen-Lake ", rmse, " ",measuredTime)

    println(τ)
end
    κ = 1.417
    θ = 0.082
    σ = 0.591
    ρ = -0.384
    v0 = 0.073
end
#
# @testset "LFShortA" begin
#     r = 0.0
#     q = 0.0
#     κ = 0.254
#     θ = 0.320
#     σ = 0.344
#     ρ = -0.557
#     v0 = 0.826
#     τ = 0.0182
#     spot = 1000.0
#
#     spot *= exp((r - q) * τ)
#     df = exp(-r * τ)
#     params = HestonParams(v0, κ, θ, ρ, σ)
#     m = 200
#     l = 8
#     refPricer = makeCosCharFuncPricer(
#         Complex,
#         Float64,
#
#         params,
#         τ,
#         2048*4,
#         32,
#     )
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64,  params, τ, m, l)
#     strike = 1400.0
#         refPrice = priceEuropean(refPricer, true, strike, spot, df)
#         price = priceEuropean(pricer, true, strike, spot, df)
#         println(
#             strike,
#             " ",
#             price,
#             " ",
#             refPrice,
#             " ",
#             price - refPrice,
#             " ",
#             price / refPrice - 1,
#         )
#
#     cinf = (params.v0+params.κ*params.θ*τ)/params.σ *sqrt(1-params.ρ^2)
#     #exp(-cinf*umax) = eps
#     #2*exp(-cinf*mmax*pi/(b-a))/(mmax*pi) = eps => exp(cinf*mmax*pi/(b-a))*mmax*pi/(b-a)*cinf = 2/eps/(b-a)*cinf => cinf*mmax*pi/(b-a) = lambertW( 2/eps/(b-a)*cinf)
#     eps = 1e-8
#     umax = -log(eps)/cinf
#     mmax = ceil(Int,umax*(pricer.b-pricer.a)/pi)+1
#     mmax = ceil(Int, sqrt(   lambertW(2*v0*τ/(eps*eps*(pricer.b-pricer.a)))/(v0*τ*pi)*(pricer.b-pricer.a)) )
#     mmax = max(32, ceil(Int, lambertW(2*cinf/(eps*(pricer.b-pricer.a)))/(cinf*pi)*(pricer.b-pricer.a)))
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64,  params, τ, mmax, l)
#             price = priceEuropean(pricer, true, strike, spot, df)
#             println(
#                 strike,
#                 " ",
#                 price,
#                 " ",
#                 refPrice,
#                 " ",
#                 price - refPrice,
#                 " ",
#                 price / refPrice - 1,
#             )
# end
#
#
# @testset "LF1Y" begin
#     r = 0.0
#     q = 0.0
#     κ = 0.1
#     θ = 0.01
#     σ = 2.0
#     ρ = 0.5
#     v0 = 0.0225
#     τ = 1.0
#     spot = 1000000.0
#
#     spot *= exp((r - q) * τ)
#     df = exp(-r * τ)
#     params = HestonParams(v0, κ, θ, ρ, σ)
#     m = 200
#     l = 10
#     refPricer = makeCosCharFuncPricer(
#         Complex,
#         Float64,
#
#         params,
#         τ,
#         2048*4,
#         32,
#     )
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64,  params, τ, m, l)
#     strike = spot*0.25
#         refPrice = priceEuropean(refPricer, false, strike, spot, df)
#         price = priceEuropean(pricer, false, strike, spot, df)
#         println(
#             strike,
#             " ",
#             price,
#             " ",
#             refPrice,
#             " ",
#             price - refPrice,
#             " ",
#             price / refPrice - 1,
#         )
#
#     cinf = (params.v0+params.κ*params.θ*τ)/params.σ *sqrt(1-params.ρ^2)
#     #exp(-cinf*umax) = eps
#     #2*exp(-cinf*mmax*pi/(b-a))/(mmax*pi) = eps => exp(cinf*mmax*pi/(b-a))*mmax*pi/(b-a)*cinf = 2/eps/(b-a)*cinf => cinf*mmax*pi/(b-a) = lambertW( 2/eps/(b-a)*cinf)
#     eps = 1e-10
#     umax = -log(eps)/cinf
#     mmax = ceil(Int,umax*(pricer.b-pricer.a)/pi)+1
#     mmax = ceil(Int, sqrt(   lambertW(2*v0*τ/(eps*eps*(pricer.b-pricer.a)))/(v0*τ*pi)*(pricer.b-pricer.a)) )
#     mmax = max(32, ceil(Int, lambertW(2*cinf/(eps*(pricer.b-pricer.a)))/(cinf*pi)*(pricer.b-pricer.a)))
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64,  params, τ, mmax, l)
#             price = priceEuropean(pricer, false, strike, spot, df)
#             println(
#                 strike,
#                 " ",
#                 price,
#                 " ",
#                 refPrice,
#                 " ",
#                 price - refPrice,
#                 " ",
#                 price / refPrice - 1,
#             )
# end
@testset "LFShort" begin
    r = 0.0
    q = 0.0
    κ = 1.0
    θ = 0.1
    σ = 1.0
    ρ = -0.9
    v0 = 0.1
    τ = 2.0 / 365
    spot = 1.0

    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    m = 256
    l = 12
    refPricer = makeCosCharFuncPricer(cf, τ, 1024, 32)
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    strikes = 1.0:0.025:1.5
    refPrices = map(strike -> priceEuropean(refPricer, true, strike, spot,  τ,df), strikes)
    for (strike,refPrice) = zip(strikes,refPrices)
        price = priceEuropean(pricer, true, strike, spot, τ, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )
        @test isapprox(refPrice, price, atol = 1e-15)
    end
    ms,tol = CharFuncPricing.findSwiftScaling(cf, τ)
    pricerS = SwiftCharFuncPricer(cf, τ, ms,l)
    for (strike,refPrice) = zip(strikes,refPrices)
        price = priceEuropean(pricerS, true, strike, spot, τ, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )        
    end
end


@testset "LFSwiftSet1" begin
    r = 0.0
    q = 0.0
    κ = 0.1
    θ = 0.01
    σ = 2.0
    ρ = 0.5
    v0 = 0.0225
    τ = 1.0
    spot = 1e6

    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    m = 256
    l = 12
    refPricer = ALCharFuncPricer(cf,  n=512)
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    strikes = [0.25e6,1e6,4e6]
    refPrices = map(strike -> priceEuropean(refPricer, strike > spot, strike, spot,  τ,df), strikes)
    for (strike,refPrice) = zip(strikes,refPrices)
        price = priceEuropean(pricer, strike > spot, strike, spot, τ, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )
        @test isapprox(refPrice, price, atol = 5e-2*refPrice)
    end
    ms,tol = CharFuncPricing.findSwiftScaling(cf, τ)
    pricerS = SwiftCharFuncPricer(cf, τ, ms,l,tol=1.0)    
    for (strike,refPrice) = zip(strikes,refPrices)
        price = priceEuropean(pricerS, strike > spot, strike, spot, τ, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )        
    end
    #vieta best.  simpson worse
    ms,tol = CharFuncPricing.findSwiftScaling(cf, τ,tol=1e-4)
    strikes = range(0.25e6,stop=4e6, length=51)
    refPrices = map(strike -> priceEuropean(refPricer, strike > spot, strike, spot,  τ,df), strikes)
    names = ["Vieta","Trapezoidal","EM1"]
    pricers = [CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,densityCalculator=CharFuncPricing.VietaDensityCalculator()),
    CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,densityCalculator=CharFuncPricing.TrapezoidalDensityCalculator()),
    CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,densityCalculator=CharFuncPricing.EulerMaclaurinDensityCalculator{1}())]
    for (name,pricer) = zip(names,pricers)
        pricesV = map(strike -> priceEuropean(pricer, strike > spot, strike, spot,  τ,df,payoffCalculator=CharFuncPricing.VietaPayoffCalculator()), strikes)
        println(name, "-Vieta ",rmsd(pricesV,refPrices)," ",maximum(abs.(pricesV-refPrices)))        
        pricesT = map(strike -> priceEuropean(pricer, strike > spot, strike, spot,  τ,df,payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{0}()), strikes)
        println(name, "-Trapezoidal ",rmsd(pricesT,refPrices)," ",maximum(abs.(pricesT-refPrices)))        
        pricesEM1 = map(strike -> priceEuropean(pricer, strike > spot, strike, spot,  τ,df,payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{1}()), strikes)
        println(name, "-EM1 ",rmsd(pricesEM1,refPrices)," ",maximum(abs.(pricesEM1-refPrices)))        
    end
end


@testset "LFSwiftSet2" begin
    r = 0.0
    q = 0.0
    κ = 1.0
    θ = 0.1
    σ = 1.0
    ρ = -0.9
    v0 = 0.1
    τ = 2.0/365
    spot = 1.0

    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    m = 256
    l = 12
    refPricer = ALCharFuncPricer(cf,  n=512)
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    strikes = [1.0,1.0064,1.064]
    refPrices = map(strike -> priceEuropean(refPricer, strike > spot, strike, spot,  τ,df), strikes)
    for (strike,refPrice) = zip(strikes,refPrices)
        price = priceEuropean(pricer, strike > spot, strike, spot, τ, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )
        @test isapprox(refPrice, price, atol = 1e-15)
    end    
    ms,tol = CharFuncPricing.findSwiftScaling(cf, τ)
    pricerS = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0)    
    for (strike,refPrice) = zip(strikes,refPrices)
        price = priceEuropean(pricerS, strike > spot, strike, spot, τ, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )        
    end
    #vieta/trap better than simpson
    ms,tol = CharFuncPricing.findSwiftScaling(cf, τ,tol=1e-4)
   strikes = range(0.9,stop=1.1, length=51)
    refPrices = map(strike -> priceEuropean(refPricer, strike > spot, strike, spot,  τ,df), strikes)
    names = ["Vieta","Trapezoidal","EM1"]
    pricers = [CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,densityCalculator=CharFuncPricing.VietaDensityCalculator()),
    CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,densityCalculator=CharFuncPricing.TrapezoidalDensityCalculator()),
    CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,densityCalculator=CharFuncPricing.EulerMaclaurinDensityCalculator{1}())]
    # pr icesV = Vector(length(names))
    for (name,pricer) = zip(names,pricers)
        pricesV = map(strike -> priceEuropean(pricer, strike > spot, strike, spot,  τ,df,payoffCalculator=CharFuncPricing.VietaPayoffCalculator()), strikes)
        println(name, "-Vieta ",rmsd(pricesV,refPrices)," ",maximum(abs.(pricesV-refPrices)))        
        pricesT = map(strike -> priceEuropean(pricer, strike > spot, strike, spot,  τ,df,payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{0}()), strikes)
        println(name, "-Trapezoidal ",rmsd(pricesT,refPrices)," ",maximum(abs.(pricesT-refPrices)))        
        pricesEM1 = map(strike -> priceEuropean(pricer, strike > spot, strike, spot,  τ,df,payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{1}()), strikes)
        println(name, "-EM1 ",rmsd(pricesEM1,refPrices)," ",maximum(abs.(pricesEM1-refPrices)))        
    end
end


#Testing https://financepress.com/2019/02/15/heston-model-reference-prices/
@testset "AlanFloat64Set" begin
    r = 0.01
    q = 0.02
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    spot = 100.0
    strikes = [80.0, 90.0, 100.0, 110.0, 120.0]
    alanPuts = [
        7.958878113256768285213263077598987193482161301733,
        12.017966707346304987709573290236471654992071308187,
        17.055270961270109413522653999411000974895436309183,
        23.017825898442800538908781834822560777763225722188,
        29.811026202682471843340682293165857439167301370697,
    ]
    alanCalls = [
        26.774758743998854221382195325726949201687074848341,
        20.933349000596710388139445766564068085476194042256,
        16.070154917028834278213466703938231827658768230714,
        12.132211516709844867860534767549426052805766831181,
        9.024913483457835636553375454092357136489051667150,
    ]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    l = 32
    m = 1024
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricer, false, strike, spot, τ, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 1e-13)
        price = priceEuropean(pricer, true, strike, spot,  τ,df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 1e-13)
    end
end
@testset "AlanJoshiFloat64Set" begin
    r = 0.01
    q = 0.02
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.01
    τ = 0.01
    spot = 100.0
    strikes = [90.0, 95.0, 100.0, 105.0, 110.0]
    alanPuts = [
        4.5183603586861772614990106188215872180542e-8,
        0.000461954855653851579672612557018857858641926937,
        0.477781171629504680023239655436072890669645669297,
        5.009501052563650299130635110520904481889436667608,
        10.008998550115123724684210555728039829315964456261,
    ]
    alanCalls = [
        9.989001595065276544935948045293485530832966049263,
        4.989963479738160122154264702582719627807098780529,
        0.467782671512844263098248405184095087949465507760,
        2.527447823194706060519991248106500619490942e-6,
        1.29932760052624920704881258510264466e-13,
    ]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    l = 32
    m = 1024
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricer, false, strike, spot, τ, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 2e-14)
        price = priceEuropean(pricer, true, strike, spot, τ,  df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 2e-14)
    end
    m,tol = CharFuncPricing.findSwiftScaling(cf, τ)
    pricerS = SwiftCharFuncPricer(cf, τ, m, 12)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricerS, false, strike, spot, τ, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 2e-7)
        price = priceEuropean(pricerS, true, strike, spot, τ,  df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 2e-7)
    end
    #Filon/Flinn require many points (Esp. Filon)
    pricerS = FlinnCharFuncPricer(cf, τ, qTol=1e-12, tTol=1e-12)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricerS, false, strike, spot, τ, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 2e-7)
        price = priceEuropean(pricerS, true, strike, spot, τ,  df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 2e-7)
    end
    
end

#=
#blackscholes like
params = HestonParams{Float64}(0.25^2, 1.0, 0.25^2, 0.0, 0.01)
cf = DefaultCharFunc(params)
refPricer = ALCharFuncPricer(cf,  n=512)
τ = 10.0
m, stol = CharFuncPricing.findSwiftScaling(cf,τ)
refPrices = map(strike -> priceEuropean(refPricer, strike > f, strike, f,  τ,df), strikes)
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j), strikes)
k=-pricerS.k1+1:pricerS.k1

narrowStrikes = strikes[end:end]
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j*128), narrowStrikes)
refUk = copy(pricerS.uk)
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j), narrowStrikes)

plot(k, refUk[ifelse.(k.>0,k,k.+2pricerS.k1)] ,label="Reference",xlab=L"k",ylab=L"V_{6,k}")
plot!(k, pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ,label="Second Euler-Maclaurin of order 0")
plot(k, pricerS.cl[ifelse.(k.>0,k,k.+2pricerS.k1)],label="",xlab=L"k", ylab=L"c_{6,k}")

prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j), narrowStrikes)
plot(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ./ refUk[ifelse.(k.>0,k,k.+2pricerS.k1)].-1),label="Second Euler-Maclaurin of order 0")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{0}()), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ./ refUk[ifelse.(k.>0,k,k.+2pricerS.k1)].-1),label="First Euler-Maclaurin of order 0")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{1}()), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ./ refUk[ifelse.(k.>0,k,k.+2pricerS.k1)].-1),label="First Euler-Maclaurin of order 1")
#prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.MidPoint)), narrowStrikes)
#plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ./ refUk[ifelse.(k.>0,k,k.+2pricerS.k1)].-1),label="Mid-point rule")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.Simpson)), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ./ refUk[ifelse.(k.>0,k,k.+2pricerS.k1)].-1),label="Simpson rule")
#prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.Boole)), narrowStrikes)
#plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] ./ refUk[ifelse.(k.>0,k,k.+2pricerS.k1)].-1),label="Boole rule")

plot!(yscale=:log10, legend=:bottomright, xlab=L"k", ylab=string("Relative error in ",L"V_{m,k}"),yticks=[1e-5,1e-3,1e-1,10])

#direct errors
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j), narrowStrikes)
plot(k, pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)],label="Second Euler-Maclaurin of order 0")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{0}()), narrowStrikes)
plot!(k, pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)],label="First Euler-Maclaurin of order 0")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{1}()), narrowStrikes)
plot!(k, pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)],label="First Euler-Maclaurin of order 1")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.MidPoint)), narrowStrikes)
plot!(k, pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)],label="Mid-point rule")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.Simpson)), narrowStrikes)
plot!(k, pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)],label="Simpson rule")
#abs errors
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j), narrowStrikes)
plot(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)]),label="Second Euler-Maclaurin of order 0")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{0}()), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)]),label="First Euler-Maclaurin of order 0")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{1}()), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)]),label="First Euler-Maclaurin of order 1")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.MidPoint)), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)]),label="Mid-point rule")
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.Simpson)), narrowStrikes)
plot!(k, abs.(pricerS.uk[ifelse.(k.>0,k,k.+2pricerS.k1)] .- refUk[ifelse.(k.>0,k,k.+2pricerS.k1)]),label="Simpson rule")

strikes=range(0.9,stop=1.1, length=101)
useVaryingJ = false
refPrices = map(strike -> priceEuropean(refPricer, strike > f, strike, f,  τ,df), strikes)
prices = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j,useVaryingJ=useVaryingJ), strikes)
prices128 = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j*128,useVaryingJ=useVaryingJ), strikes)
pricesEM0 = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j,useVaryingJ=useVaryingJ, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{0}()),strikes)
pricesEM1 = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j,useVaryingJ=useVaryingJ, payoffCalculator=CharFuncPricing.EulerMaclaurinPayoffCalculator{1}()), strikes)
pricesS = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j,useVaryingJ=useVaryingJ, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.Simpson)), strikes)


plot(strikes, abs.(prices-refPrices),label="Second Euler-Maclaurin of order 0" )
plot!(strikes, abs.(pricesEM0-refPrices),label="First Euler-Maclaurin of order 0" )
plot!(strikes, abs.(pricesEM1-refPrices),label="First Euler-Maclaurin of order 1" )
plot!(strikes, abs.(prices128-refPrices),label="Large J limit" )
plot!(strikes, abs.(pricesS-refPrices),label="Simpson rule" )
plot!(xlab="Strike",ylab="Error in option price",yscale=:log10,legend=:bottomleft)

using LaTeXStrings
ms=6;l=8;
pricerV = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0, densityCalculator=CharFuncPricing.VietaDensityCalculator())
pricerV8 = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,l,tol=1.0,Jfactor=8)
pricerT = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms, l, tol=1.0, densityCalculator=CharFuncPricing.TrapezoidalDensityCalculator())
pricerEM1 = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms, l, tol=1.0, densityCalculator=CharFuncPricing.EulerMaclaurinSlowDensityCalculator{1}())
 k=-pricerV8.k1+1:pricerV8.k1
 plot(k,pricerV.cl[ifelse.(k.>0,k,k.+2pricerV8.k1)]-pricerV8.cl[ifelse.(k.>0,k,k.+2pricerV8.k1)],label="Vieta")
plot!(k,pricerT.cl[ifelse.(k.>0,k,k.+2pricerV8.k1)]-pricerV8.cl[ifelse.(k.>0,k,k.+2pricerV8.k1)],label="Trapezoidal")
plot!(k,pricerEM1.cl[ifelse.(k.>0,k,k.+2pricerV8.k1)]-pricerV8.cl[ifelse.(k.>0,k,k.+2pricerV8.k1)],label="Euler Maclaurin order 1")
plot!(ylab=string("Error in ",L"c_{6,k}"),xlab=L"k")
plot!(size=(400,400),margins=3Plots.mm)

refPrices = map(strike -> priceEuropean(refPricer, strike > f, strike, f,  τ,df), strikes)
prices = map(strike -> priceEuropean(pricerV, strike > f, strike, f,  τ,df,J=pricerV.j*16), strikes)
prices8 = map(strike -> priceEuropean(pricerV8, strike > f, strike, f,  τ,df,J=pricerV.j*16), strikes)
pricesT = map(strike -> priceEuropean(pricerT, strike > f, strike, f,  τ,df,J=pricerV.j*16), strikes)
pricesEM1 = map(strike -> priceEuropean(pricerEM1, strike > f, strike, f,  τ,df,J=pricerV.j*16), strikes)
#pricesS = map(strike -> priceEuropean(pricerS, strike > f, strike, f,  τ,df,J=pricerS.j,useVaryingJ=useVaryingJ, payoffCalculator=CharFuncPricing.IntegralPayoffCalculator(CharFuncPricing.Simpson)), strikes)

plot(strikes, abs.(prices-refPrices),label="Second Euler-Maclaurin of order 0" )
plot!(strikes, abs.(pricesT-refPrices),label="First Euler-Maclaurin of order 0" )
plot!(strikes, abs.(pricesEM1-refPrices),label="First Euler-Maclaurin of order 1" )
plot!(strikes, abs.(prices8-refPrices),label="Large J limit" )
#plot!(strikes, abs.(pricesS-refPrices),label="Simpson rule" )
plot!(xlab="Strike",ylab="Error in option price",legend=:bottomright)


=#