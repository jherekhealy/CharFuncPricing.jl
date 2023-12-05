using CharFuncPricing, Test
@testset "BroadieKayaSVCJ" begin
    hestonParams = HestonParams(0.007569, 3.46, 0.008, -0.82, 0.14)
    params = SVCJParams(hestonParams, 0.47, -0.1, 0.0001^2, 0.05, -0.38)
    df = exp(-0.0319)
    cf = DefaultCharFunc(params)
    refPrice = 6.8619 # from Broadie Kaya paper.
    pricer160 = CharFuncPricing.AdaptiveFilonCharFuncPricer(cf,1.0,qTol=1e-8,myTrans=CharFuncPricing.IdentityTransformation(0.0,160.0))
    price = priceEuropean(pricer160,true, 100.0, 100.0/df,1.0,df)
    @test isapprox(refPrice,price, atol=1e-4)
    pricerCos = CharFuncPricing.makeCosCharFuncPricer(cf,1.0,1024,16)
    price = priceEuropean(pricerCos,true, 100.0, 100.0/df,1.0,df)
    @test isapprox(refPrice,price, atol=1e-4)
    p = JoshiYangCharFuncPricer(cf,1.0)
    price = priceEuropean(p,true, 100.0, 100.0/df,1.0,df)
    @test isapprox(refPrice,price, atol=1e-4)
    
    #adaptive filon is poor for large strikes in general or very small strikes with AL transf. Cos stays accurate on this example
end