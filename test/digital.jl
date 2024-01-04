using CharFuncPricing, Test

@testset "HestonDigital" begin
    hestonParams = HestonParams(0.007569, 3.46, 0.008, -0.82, 0.14)
    cf = DefaultCharFunc(hestonParams)
    τ = 1.0
    pricerCos = CharFuncPricing.makeCosCharFuncPricer(cf, τ,1024,16)
    pal = ALCharFuncPricer(cf)
    ms, tol = CharFuncPricing.findSwiftScaling(cf, τ)
    pricerS = CharFuncPricing.makeSwiftCharFuncPricer(cf, τ, ms,16)
    strikes = [0.7,1.0,1.3]
    forward = 1.05
    for strike = strikes
        p = priceDigital(pricerCos,false,strike,forward,τ,1.0)
        pd = (priceEuropean(pal,false, strike+1e-4, forward, τ,1.0)-priceEuropean(pal,false, strike-1e-4, forward, τ,1.0))/2e-4
        ps = priceDigital(pricerS, false, strike, forward, τ, 1.0)
        println(strike," ", p, " ", pd-p, " ",ps, " ",ps-p)
        @test isapprox(p, pd, atol=3e-7)
        @test isapprox(ps, p, atol=1e-8)
        p = priceDigital(pricerCos,true,strike,forward,τ,1.0)
        pd = -(priceEuropean(pal,true, strike+1e-4, forward, τ,1.0)-priceEuropean(pal, true, strike-1e-4, forward, τ,1.0))/2e-4
        ps = priceDigital(pricerS, true, strike, forward, τ, 1.0)
        @test isapprox(p, pd, atol=3e-7)
        @test isapprox(ps, p, atol=1e-8)
    end

end
 
