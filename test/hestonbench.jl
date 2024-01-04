using CharFuncPricing, Test, Printf
using StatsBase

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
    isCall = true
    paramset = Vector(undef, 0)
    for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
    params = HestonParams(v0, κ, θ, ρ, σ)
    for strike in strikes
        push!(paramset, (τ, strike, params))
    end
end

    pricesA1 = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward, τ, 1.0)
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
            pricer = ALCharFuncPricer(cf, n=1000)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward, τ, 1.0)
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
            pricer = ALCharFuncPricer(cf,  n=2000)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward, τ, 1.0)
                pricesA3[counter] = price
            end
        end
    end

    pricesA0GL = zeros(length(paramset));
    qgl = ModlobQuadrature(1e-10)
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, qgl)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward,τ, 1.0)
                pricesA0GL[counter] = price
            end
        end
    end

    pricesA0L = zeros(length(paramset));
    qde = DEQuadrature(1e-10)
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, qde)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward,τ, 1.0)
                pricesA0L[counter] = price
            end
        end
    end
    pricesS1 = Vector{Float64}(undef, 0);
    elapsed = @elapsed begin
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
        params = HestonParams(v0, κ, θ, ρ, σ)
        cf = DefaultCharFunc(params)
        m, tol = CharFuncPricing.findSwiftScaling(cf, τ, tol=1e-4)
        refPricer = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, m, 2, tol=max(2tol,1e-4))
        for strike in strikes
            price = priceEuropean(refPricer, isCall, strike, forward, τ, 1.0)
            push!(pricesS1, price)
        end
    end
    end

    pricesJY1 =  Vector{Float64}(undef, 0);
    elapsed = @elapsed begin
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
        params = HestonParams(v0, κ, θ, ρ, σ)
        cf = DefaultCharFunc(params)
        refPricer = CharFuncPricing.JoshiYangCharFuncPricer(cf, τ)
        for strike in strikes
            price = priceEuropean(refPricer, isCall, strike, forward, τ, 1.0)
            push!(pricesJY1, price)
        end
    end
    end

    pricesC1 = Vector{Float64}(undef, 0);
    elapsed = @elapsed begin
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
        params = HestonParams(v0, κ, θ, ρ, σ)
        cf = DefaultCharFunc(params)
        refPricer = makeCosCharFuncPricer(cf, τ, 200 , 8)
        for strike in strikes
            price = priceEuropean(refPricer, isCall, strike, forward, τ, 1.0)
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
            price = priceEuropean(refPricer, isCall, strike, forward, τ, 1.0)
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
                price = priceEuropean(pricer, isCall, strike, forward, 1.0)
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
                price = priceEuropean(pricer, isCall, strike, forward, 1.0)
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
                price = priceEuropean(pricer, isCall, strike, forward, 1.0)
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
            pricer = FlinnCharFuncPricer(cf, τ, tTol=1e-6, qTol=1e-10)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward, τ, 1.0)
                pricesF2B[counter] = price
            end
        end
    end

    pricesJ = zeros(length(paramset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = JoshiYangCharFuncPricer(cf, τ, n=64)
            for strike in strikes
                counter += 1
                price = priceEuropean(pricer, isCall, strike, forward, τ, 1.0)
                pricesJ[counter] = price
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
end
#
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


    quadrature = TanhSinhQuadrature(200,eps())
    fpricesA1 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf,quadrature)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, τ, 1.0)
                fpricesA1[counter] = price
            end
        end
    end

    quadrature = TanhSinhQuadrature(1000,eps())
fpricesA2 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, quadrature)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, τ, 1.0)
                fpricesA2[counter] = price
            end
        end
    end

    quadrature = TanhSinhQuadrature(2000,eps())
    fpricesA3 = zeros(length(fparamset));
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, quadrature)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, τ, 1.0)
                fpricesA3[counter] = price
            end
        end
    end

    quadrature = ModlobQuadrature(1e-10)
    fpricesA0L = zeros(length(fparamset))
    elapsed = @elapsed begin
        local counter = 0
        for τ in τs, v0 in v0s, θ in θs, κ in κs, σ in σs, ρ in ρs
            params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
            cf = DefaultCharFunc(params)
            pricer = ALCharFuncPricer(cf, quadrature)
            for forward in forwards
                counter += 1
                price = priceEuropean(pricer, false, strike, forward, τ, 1.0)
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
end
