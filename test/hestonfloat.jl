using CharFuncPricing, Test
using TaylorSeries
using StatsBase

@testset "HestonLong" begin
τ=10.0
forward= 10000.0
params=HestonParams(1.0, 2.0, 0.0025, 0.5, 0.1)
strikes = [100.0001, 101.0, 110.0, 200.0, 1000.0, 10000.0]
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, τ, n=0000), makeCosCharFuncPricer(cf, τ, 200, 8), makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10), FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, τ, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end



τ=10.0
forward = 100.0
params = HestonParams(0.0001, 0.1, 0.25, 0.95, 3.0)
strike = 100.0001
# 0.16558202265851305 13.845748912814255
# 10.0 0.0001 0.1 0.25 0.95 3.0 101.0 1.138393784200611 14.617071148832778
# 10.0 0.0001 0.1 0.25 0.95 3.0 110.0 10.036664633388256 22.41077792963597
# 10.0 0.0001 0.1 0.25 0.95 3.0 200.0 100.00004778409894 112.733332472877
# 10.0 0.0001 0.1 0.25 0.95 3.0 1000.0 899.9999999999567 911.8681921880134
# 10.0 0.0001 0.1 0.25 0.95 3.0 10000.0 9900.0 9917.669236275275
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 8)
cosPricer2 = makeCosCharFuncPricer(cf, τ, 200, 8)
cosPricer10 = makeCosCharFuncPricer(cf, τ, 1000, 8)
fPricer = FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=200)
isCall = false
#bad price for adapt-cos vs cos 8
for strike in strikes
    count += 1
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceC2 = priceEuropean(cosPricer2, isCall, strike, forward, 1.0)
    priceC10 = priceEuropean(cosPricer10, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," & ",priceC," & ",priceC2," & ",priceC10," & ",priceF, " & ",priceA," \\\\")
end
#10.0 0.0001 0.01 1.0 -0.95 3.0 10000.0 9900.0 9882.675792888303
params = HestonParams(0.0001, 0.01, 1.0, -0.95, 3.0)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 24)
fPricer = FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=200)
isCall = false
for strike in strikes
    count += 1
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end

τ=10.0
strike=101.0
params=HestonParams{Float64}(0.0001, 0.01, 0.0001, 0.95, 3.0)
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, τ, n=0000), makeCosCharFuncPricer(cf, τ, 200, 8), makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10), FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, τ, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end


τ=10.0
strike=101.0
params= HestonParams{Float64}(1.0, 0.01, 0.25, 0.95, 3.0)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 8)
fPricer = FlinnCharFuncPricer(HestonCVCharFunc(cf), τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=200)
isCall = false
for strike in strikes
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end

#(30.0, 100.0001, HestonParams{Float64}(1.0, 0.1, 1.0, 0.95, 0.1))
τ=30.0
params =  HestonParams{Float64}(1.0, 0.1, 1.0, 0.95, 0.1)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 24)
fPricer = FlinnCharFuncPricer(HestonCVCharFunc(cf), τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=200)
isCall = false
for strike in strikes
    count += 1
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end


τ=0.0025
params = HestonParams{Float64}(0.0001, 0.01, 0.0001, -0.95, 3.0)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 8)
fPricer = FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=10000)
isCall = false
for strike in strikes
    count += 1
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end

#flinn breaks strongly with or without cv. Also due to bad truncatino range (better if we use short mat criteria as well)
τ=0.5
strike = 100.0
params = HestonParams{Float64}(0.0025, 2.0, 0.25, 0.5, 0.0001)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 8)
fPricer = FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=10000)
isCall = false
for forward in strikes
    count += 1
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end

#DE adaptive largest error abs
τ= 0.1
strike= 101.0
forward=100.0
params = HestonParams{Float64}(0.0001, 2.0, 0.04, -0.5, 3.0)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 8)
fPricer = FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10)
aPricer = ALCharFuncPricer(cf, τ, n=10000)
isCall = false
for strike in strikes
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end

#TFlinn largest error with 1e-8 (ok with 1e-10)
params = paramset[2497][3]
HestonParams{Float64}(0.0001, 2.0, 0.04, 0.0, 3.0)

julia> τ=paramset[2497][1]
0.0025

#Largest error of TFlinn 1e-10 and AL-200. Cos requires million points!
#Interestingly, AL-10000 is likley wrong as AL-0 with tol=1e-20 = TFLinn price, strike = 101 & 200.
τ=30.0
forward = 100.0
params = HestonParams{Float64}(0.0025, 0.1, 0.0001, 0.1, 3.0)
cf = DefaultCharFunc(params)
cosPricer = makeCosCharFuncPricer(cf, τ, 16)
fPricer = FlinnCharFuncPricer(cf, τ, qTol=1e-15)
aPricer = ALCharFuncPricer(cf, τ, n=10000)
isCall = false
for strike in strikes
    priceC = priceEuropean(cosPricer, isCall, strike, forward, 1.0)
    priceF = priceEuropean(fPricer, isCall, strike, forward, 1.0)
    priceA = priceEuropean(aPricer, isCall, strike, forward, 1.0)
    println(strike," ",priceC," ",priceF," ",priceA)
end
#Largest error on AL-1000 vs 2000
τ=10.0
forward = 100.0
params = HestonParams{Float64}(0.0001, 0.01, 0.0001, 0.95, 1.0)
#(10.0, 100.0001, HestonParams{Float64}(0.0001, 0.01, 0.0001, 0.95, 1.0))
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, τ, n=0000), makeCosCharFuncPricer(cf, τ, 200, 8), makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10), FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, τ, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, x, forward, 1.0),strikes)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, x, forward, 1.0), strikes)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end
#largest error n=200
τ=2.0
strike = 100.0
forward = 200.0
params = HestonParams{Float64}(0.0001, 2.0, 0.0001, 0.95, 1.0)
cf = DefaultCharFunc(params)
pricers = [ALCharFuncPricer(cf, τ, n=0000), makeCosCharFuncPricer(cf, τ, 200, 8), makeCosCharFuncPricer(cf, τ, 8, tol=1e-8), FlinnCharFuncPricer(cf, τ, tTol = 1e-10, qTol=1e-10), FlinnCharFuncPricer(cf, τ, tTol = 1e-40, qTol=1e-10), AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-8), ALCharFuncPricer(cf, τ, n=200)]
pricerNames = ["Reference", "Cos (M=200)", "Cos-Adaptive", "Flinn-Truncated (1e-8)", "Flinn-Truncated (1e-40)", "Flinn-Transformed", "Andersen-Lake (n=200)"]
isCall = false
#bad price for adapt-cos vs cos 8
refPrices = map(x-> priceEuropean(pricers[1], isCall, strike, x, 1.0),forwards)
for price in refPrices
    @printf("& %.8e", price)
end
println()
for (pricer,name) in zip(pricers,pricerNames)
    prices = map(x -> priceEuropean(pricer, isCall, strike, x, 1.0), forwards)
    print(name," ")
    for (price,refPrice) in zip(prices,refPrices)
        @printf("& %.2e", price - refPrice)
    end
    println("\\\\")
    # println(strike," & ",priceC," & ",priceC2," & ",priceF," & ",priceF2, " & ",priceA," \\\\")
end

end

# bench test from AndersenLake (goal: verify accuracy of Cos & Flinn)
@testset "HestonALBench" begin
    forward = 100.0
    strikes = [100.0001, 101.0, 110.0, 200.0, 1000.0, 10000.0]
    τs = [0.0025, 0.1, 0.5, 2.0, 10.0, 30.0]
    v0s = [0.0001, 0.0025, 0.04, 0.25, 1.0]
    θs = [0.0001, 0.0025, 0.04, 0.25, 1.0]
    κs = [0.01, 0.1, 0.5, 2.0]
    σs = [0.0001, 0.1, 0.5, 1.0, 3.0]
    ρs = [-0.95, -0.5, -0.1, 0.0, 0.1, 0.5, 0.95]
    paramset = Vector(undef, 0)
    for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
    params = HestonParams(v0, κ, θ, ρ, σ)
    for strike in strikes
        push!(paramset, (τ, strike, params))
    end
end

    pricesC1 = Vector{Float64}(undef, 0);
    elapsed = @elapsed begin
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
        params = HestonParams(v0, κ, θ, ρ, σ)
        cf = DefaultCharFunc(params)
        refPricer = makeCosCharFuncPricer(cf, τ, 200 , 8)
        for strike in strikes
            price = priceEuropean(refPricer, false, strike, forward, 1.0)
            push!(pricesC1, price)
        end
    end
    end

    pricesC2 = Vector{Float64}(undef, 0);
    elapsed = @elapsed begin
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
        params = HestonParams(v0, κ, θ, ρ, σ)
        cf = DefaultCharFunc(params)
        refPricer = makeCosCharFuncPricer(cf, τ, 8)
        for strike in strikes
            price = priceEuropean(refPricer, false, strike, forward, 1.0)
            push!(pricesC2, price)
        end
    end
end

    pricesF1 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = CharFuncPricing.AdaptiveFlinnCharFuncPricer(cf, τ)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesF1[counter] = price
            end
        end
    end

    pricesF2 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = CharFuncPricing.AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-9)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesF2[counter] = price
            end
        end
    end


    pricesF3 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = CharFuncPricing.AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-10)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesF3[counter] = price
            end
        end
    end

    pricesF2B = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = FlinnCharFuncPricer(cf, τ, tTol=1e-10, qTol=1e-10)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesF2B[counter] = price
            end
        end
    end


    pricesA1 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesA1[counter] = price
            end
        end
    end

    pricesA2 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ, n=1000)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesA2[counter] = price
            end
        end
    end

    pricesA3 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ, n=2000)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesA3[counter] = price
            end
        end
    end

    pricesA0L = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ, n=0)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                pricesA0L[counter] = price
            end
        end
    end

    allp = [pricesC1, pricesC2, pricesF1, pricesF2, pricesF3, pricesF2B, pricesA1, pricesA2, pricesA0L]
    names = ["Cos","Cos-Adaptive","Flinn-Transf.","Flinn-Transf. (1e-9)", "Flinn-Transf. (1e-10)","Flinn (1e-10)","Andersen-Lake", "Andersen-Lake 1000", "Andersen-Lake GL"]
    for (p,name) in zip(allp,names)
        @printf("%s & %2e & %2e & %2e & %2e\n",name, maxad(p, pricesA3), rmsd(p, pricesA3),  maxad(p ./ pricesA3, ones(length(pricesA3))), rmsd(p ./ pricesA3, ones(length(pricesA3))))
    end

    allp = [pricesA1, pricesA2, pricesA0L]
    names = ["Andersen-Lake", "Andersen-Lake 1000", "Andersen-Lake GL"]
    for (p,name) in zip(allp,names)
        p3 = pricesA3
        @printf("%s & %.1e & %.1e & %.1e & %.1e\n",name, maxad(p, p3), rmsd(p, p3),  maxad(p ./ pricesA3, ones(length(p3))), rmsd(p ./ p3, ones(length(p3))))
    end
    rmsd(pricesC1, pricesA3)
    mre = maxad(pricesA1 ./ pricesA2, ones(length(pricesA1))
    rrmse = rmsd(pricesA1 ./ pricesA2, ones(length(pricesA1)))

    thrIndices = findall(z -> z == 0, pricesA0)
    thrIndices100 = findall(z -> z > 1e-100, pricesA3)
    thrIndices = findall(z -> z > 1e-8, pricesA0)
    armse = rmsd(pricesA1, pricesA2)
    mre = maxad(pricesA1 ./ pricesA2, ones(length(pricesA1))
    rrmse = rmsd(pricesA1 ./ pricesA2, ones(length(pricesA1)))

end

@testset "HestonALBenchForward" begin
    strike = 100.0
    forwards = [100.0,100.0001, 101.0, 110.0, 200.0, 1000.0, 10000.0]
    τs = [0.0025, 0.1, 0.5, 2.0, 10.0, 30.0]
    v0s = [0.0001, 0.0025, 0.04, 0.25, 1.0]
    θs = [0.0001, 0.0025, 0.04, 0.25, 1.0]
    κs = [0.01, 0.1, 0.5, 2.0]
    σs = [0.0001, 0.1, 0.5, 1.0, 3.0]
    ρs = [-0.95, -0.5, -0.1, 0.0, 0.1, 0.5, 0.95]
    fparamset = Vector(undef, 0)
    for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
    params = HestonParams(v0, κ, θ, ρ, σ)
    for forward in forwards
        push!(fparamset, (τ, forward, params))
    end
end

fpricesC1 = Vector{Float64}(undef, 0);
elapsed = @elapsed begin
    for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
    params = HestonParams(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    refPricer = makeCosCharFuncPricer(cf, τ, 200 , 8)
    for forward in forwards
        price = priceEuropean(refPricer, false, strike, forward, 1.0)
        push!(fpricesC1, price)
    end
end
end

fpricesC2 = Vector{Float64}(undef, 0);
elapsed = @elapsed begin
    for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
    params = HestonParams(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    refPricer = makeCosCharFuncPricer(cf, τ, 8)
    for forward in forwards
        price = priceEuropean(refPricer, false, strike, forward, 1.0)
        push!(fpricesC2, price)
    end
end
end
    fpricesF1 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = AdaptiveFlinnCharFuncPricer(cf, τ)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesF1[counter] = price
            end
        end
    end

    fpricesF2 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = CharFuncPricing.AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-9)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesF2[counter] = price
            end
        end
    end

    fpricesF3 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = CharFuncPricing.AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-10)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesF3[counter] = price
            end
        end
    end
    fpricesF2B = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = FlinnCharFuncPricer(cf, τ, tTol=1e-10, qTol=1e-10)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesF2B[counter] = price
            end
        end
    end


    fpricesA1 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesA1[counter] = price
            end
        end
    end

    fpricesA2 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ, n=1000)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesA2[counter] = price
            end
        end
    end

    fpricesA3 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ, n=2000)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesA3[counter] = price
            end
        end
    end

    fpricesA0L = zeros(length(fparamset))
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, τ, n=0)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, 1.0)
                fpricesA0L[counter] = price
            end
        end
    end

    allp = [fpricesC1, fpricesC2, fpricesF1, fpricesF2, fpricesF2B, fpricesA1, fpricesA2, fpricesA0L]
    names = ["Cos","Cos-Adaptive","Flinn-Transf.","Flinn-Transf. (1e-10)","Flinn (1e-10)","Andersen-Lake", "Andersen-Lake 1000", "Andersen-Lake GL"]
    for (p,name) in zip(allpf,names)
        @printf("%s & %2e & %2e & %2e & %2e\n",name, maxad(p, fpricesA3), rmsd(p, fpricesA3),  maxad(p ./ fpricesA3, ones(length(fpricesA3))), rmsd(p ./ fpricesA3, ones(length(fpricesA3))))
    end

    allp = [pricesC1, pricesC2, pricesF1, pricesF3, pricesF2B, pricesA1, pricesA2, pricesA0L]
    allpf = [fpricesC1, fpricesC2, fpricesF1, fpricesF3, fpricesF2B, fpricesA1, fpricesA2, fpricesA0L]
    names = ["Cos","Cos-Adaptive","Flinn-Transf.","Flinn-Transf. (1e-10)","Flinn (1e-10)","Andersen-Lake", "Andersen-Lake 1000", "Andersen-Lake GL"]
    for (pk, pf, name) in zip(allp, allpf, names)
        p3 = vcat(pricesA3,fpricesA3)
        ind = findall(z -> z > 1e-90, p3)
        p = vcat(pk, pf)
        @printf("%s & %.1e & %.1e & %.1e & %.1e\n",name, maxad(p[ind], p3[ind]), rmsd(p[ind], p3[ind]),  maxad(p[ind] ./ p3[ind], ones(length(p3[ind]))), rmsd(p[ind] ./ p3[ind], ones(length(p3[ind]))))
    end

    allpf = [fpricesA1, fpricesA2]
    names = ["Andersen-Lake", "Andersen-Lake 1000"]
    for (pf, name) in zip(allpf, names)
        p3 = fpricesA3
        ind = findall(z -> z > 1e-90, p3)
        p =  pf
        @printf("%s & %.1e & %.1e & %.1e & %.1e\n",name, maxad(p[ind], p3[ind]), rmsd(p[ind], p3[ind]),  maxad(p[ind] ./ p3[ind], ones(length(p3[ind]))), rmsd(p[ind] ./ p3[ind], ones(length(p3[ind]))))
    end


    thrIndices = findall(z -> z == 0, pricesA0)
    armse = rmsd(pricesA1, pricesA2)
    mre = maxad(pricesA1 ./ pricesA2, ones(length(pricesA1))
    rrmse = rmsd(pricesA1 ./ pricesA2, ones(length(pricesA1)))

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
    refPrice = priceEuropean(refPricer, false, strike, forward, 1.0)
    price = priceEuropean(pricer, false, strike, forward, 1.0)
    @test isapprox(refPrice, price, atol = 5e-4)
    @test isapprox(209.820637, price, atol = 1e-5)

    ccf = HestonCVCharFunc(cf)
    pricer = FlinnCharFuncPricer(ccf, τ, qTol = tol, tTol = tol)
    price = priceEuropean(pricer, false, strike, forward, 1.0)
    @test isapprox(refPrice, price, atol = 1e-4)
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
        refPrice = priceEuropean(refPricer, false, strike, spot, df)
        price = priceEuropean(pricer, false, strike, spot, df)
        # println(
        #     strike,
        #     " ",
        #     price,
        #     " ",
        #     refPrice,
        #     " ",
        #     price - refPrice,
        #     " ",
        #     price / refPrice - 1,
        # )
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
        price = priceEuropean(pricer, false, strike, spot, df)
        refPrices[i] = price
    end
    end
    println("Ref ", measuredTime)

    measuredTime = @benchmark begin
        pricer = makeCosCharFuncPricer(cf, τ, 200,8)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, df)
        prices[i] = price
    end
    end
    rmse = rmsd(prices, refPrices)
    println("Cos ", rmse, " ",measuredTime)

    measuredTime = @benchmark begin
        pricer = FlinnCharFuncPricer(cf, τ)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, df)
        prices[i] = price
    end
    end
    rmse = rmsd(prices, refPrices)
    println("Truncated Flinn ", rmse, " ",measuredTime)

    measuredTime = @benchmark begin
        pricer = AdaptiveFlinnCharFuncPricer(cf, τ)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, df)
        prices[i] = price
    end
    end
    rmse = rmsd(prices, refPrices)
    println("Transformed Flinn ", rmse, " ",measuredTime)

    measuredTime = @benchmark begin
        pricer = ALCharFuncPricer(cf, τ,n=200)
        for (i, strike) in enumerate(0.4:0.012:1.6)
        price = priceEuropean(pricer, false, strike, spot, df)
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
    for strike = 1.0:0.025:1.5
        refPrice = priceEuropean(refPricer, true, strike, spot, df)
        price = priceEuropean(pricer, true, strike, spot, df)
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
        price = priceEuropean(pricer, false, strike, spot, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 1e-13)
        price = priceEuropean(pricer, true, strike, spot, df)
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
        price = priceEuropean(pricer, false, strike, spot, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 2e-14)
        price = priceEuropean(pricer, true, strike, spot, df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 2e-14)
    end
end
@testset "CumulantsJoshi" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.01
    τ = 0.01
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{Float64},Taylor1{Complex},Type}(params, Complex)
    t = Taylor1(Float64, 4)
    cft = CharFuncPricing.evaluateLogCharFuncCui(cf, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-16)
end
@testset "CumulantsAlan" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{Float64},Taylor1{Complex},Type}(params, Complex)
    t = Taylor1(Float64, 4)
    cft = CharFuncPricing.evaluateLogCharFuncCui(cf, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-15)
end
@testset "CumulantsNearZeroKappa" begin
    #Low accuracy with small kappa due to powers of kappa in the denominator.
    κ = 0.02
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{Float64},Taylor1{Complex},Type}(params, Complex)
    t = Taylor1(Float64, 4)
    cft = CharFuncPricing.evaluateLogCharFuncCui(cf, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-3)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-3)
end
