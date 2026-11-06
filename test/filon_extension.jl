using Test, CharFuncPricing

#regressions of the extension of the range of the Filon and Flinn quadratures

@testset "modsim-exact-endpoints" begin
    b0 = 3.3
    a, b, c = 2b0, 4b0, 8b0
    @test (a + b) / 2 - (b - a) / 2 != a #the end point is not recovered from the center and half width
    points = Float64[]
    f(x) = (push!(points, x); cos(x) * exp(-x / 5))
    CharFuncPricing.modsim(f, a, b, 1e-8)
    @test minimum(points) == a
    @test maximum(points) == b
    CharFuncPricing.modsim(f, b, c, 1e-8)
    #the point shared by the adjacent intervals must be a single key of the grid
    @test count(x -> abs(x - b) < 1e-10, unique(points)) == 1
end

@testset "modsim-integralEstimate" begin
    n = Ref(0)
    f(x) = (n[] += 1; 1e-10 * sin(40x))
    exact = 1e-10 * (1 - cos(120.0)) / 40
    n[] = 0
    relative = CharFuncPricing.modsim(f, 0.0, 3.0, 1e-8)
    nRelative = n[]
    n[] = 0
    estimated = CharFuncPricing.modsim(f, 0.0, 3.0, 1e-8, integralEstimate=1.0)
    #the tolerance is relative to the integral so far and not to the small integral on the new interval
    @test n[] < nRelative
    @test isapprox(relative, exact, rtol=1e-6)
    @test isapprox(estimated, exact, atol=1e-8)
end

@testset "Filon-Flinn-extension" begin
    #(params, τ, forward): the truncation rule alone is too short, too long for the first panels to see the integrand
    cases = [(HestonParams(1.0, 2.0, 0.0025, 0.5, 0.1), 10.0, 10000.0),
        (HestonParams(1.0, 2.0, 1e-4, 0.0, 1e-4), 30.0, 100.0),
        (HestonParams(0.0001, 2.0, 0.04, -0.5, 3.0), 0.1, 100.0)]
    for (params, τ, forward) in cases
        cf = DefaultCharFunc(params)
        reference = ALCharFuncPricer(cf, n=2000)
        strikes = (forward / 2, forward, 1.1 * forward)
        error(p) = maximum(abs(priceEuropean(p, true, k, forward, τ, 1.0) - priceEuropean(reference, true, k, forward, τ, 1.0)) for k in strikes)
        for pricer in (FlinnCharFuncPricer, FilonCharFuncPricer)
            short = pricer(cf, τ, tTol=1e-10, qTol=1e-10, b=2.0, maxTailIterations=0)
            @test short.b == 2.0
            @test error(short) > 1e-2
            @test pricer(cf, τ, tTol=1e-10, qTol=1e-10, b=2.0, maxTailIterations=3).b == 8.0 #at most maxTailIterations - 1 doublings
            @test error(pricer(cf, τ, tTol=1e-10, qTol=1e-10)) < 1e-5
        end
        #the Filon rule pairs the sorted points by three: an end point stored twice shifts the pairs on the rest of the range
        points = FilonCharFuncPricer(cf, τ, tTol=1e-10, qTol=1e-10, b=3.3).kcos[1, :]
        @test isodd(length(points))
        @test issorted(points)
        @test minimum(diff(points)) > 1e-10
    end
end

@testset "Filon-extension-not-finite" begin
    #the Adams scheme with few steps is not finite for large u: the extension stops at the last finite point
    params = RoughHestonParams(0.05, -0.65, 0.4089, t -> 0.025)
    cf = AdamsCharFunc(params, 100)
    τ = 1.0
    pricer = FilonCharFuncPricer(cf, τ, tTol=1e-8, qTol=1e-6, b=10.0)
    @test isfinite(pricer.b)
    @test pricer.b < 100
    @test isapprox(priceEuropean(pricer, true, 1.0, 1.0, τ, 1.0), priceEuropean(FCCCharFuncPricer(cf, τ, qTol=1e-8), true, 1.0, 1.0, τ, 1.0), atol=1e-3)
end

@testset "AdaptiveFlinn-minDepth" begin
    #with ρ = 0 the sine integrand vanishes on the first panels, which pass the error test without resolving the cosine one
    cf = DefaultCharFunc(HestonParams(1e-4, 2.0, 0.25, 0.0, 1.0))
    τ = 2.0
    reference = priceEuropean(ALCharFuncPricer(cf, n=2000), true, 100.0, 100.0, τ, 1.0)
    @test isapprox(priceEuropean(AdaptiveFlinnCharFuncPricer(cf, τ), true, 100.0, 100.0, τ, 1.0), reference, atol=1e-7)
end
