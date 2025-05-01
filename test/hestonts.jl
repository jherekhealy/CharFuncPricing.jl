using CharFuncPricing
@testset "HestonTS-BGM-Table3.16" begin
    Tstart = collect(0:39) ./ 4;
    theta = collect( 0.04 .+ (0:39)*0.05/100);
    sigma = collect(0.30 .+ (0:39)*0.5/100);
    rho = collect(-0.2 .+ (0:39)*0.35/100);
    kappa = ones(40) .* 3.0;
    v0 = 0.04
    spot = 100.0
    strike = 100.0
    r = 0.0
    q = 0.0
    
    refPrices100 = [3.93, 5.53, 7.85, 11.23, 13.92, 18.37, 22.15, 27.17]
    params = HestonTSParams(v0, kappa, theta, rho, sigma, Tstart)
    cf = DefaultCharFunc(params)
    maturities = [3/12, 6/12, 1, 2, 3, 5, 7, 10];
    prices = zeros(length(maturities))
    pricesCos = zeros(length(maturities))
    for (i,τ) = enumerate(maturities)
        df = exp(-r*τ)
        pricer = JoshiYangCharFuncPricer(cf, τ)
        prices[i] = priceEuropean(pricer, false, strike, spot*exp((r-q)*τ), τ, df)
        pricerCos = makeCosCharFuncPricer(cf, τ, 128, 8)
        pricesCos[i] = priceEuropean(pricerCos, false, strike, spot, τ, df)

    end
    @test isapprox(refPrices100,prices, atol=1e-2)
    @test isapprox(refPrices100,pricesCos, atol=1e-2)

end