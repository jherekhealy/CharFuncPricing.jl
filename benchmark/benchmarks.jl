#Benchmark suite, in the PkgBenchmark layout (defines SUITE). Run with
#   julia --project=benchmark benchmark/run.jl
using CharFuncPricing, BenchmarkTools

const SUITE = BenchmarkGroup()

let
    τ = 1.0
    forward = 100.0
    cf = DefaultCharFunc(HestonParams(0.04, 4.0, 0.25, -0.5, 1.0))
    strikes = collect(range(50.0, 200.0, length=100))
    isCalls = falses(length(strikes))
    priceAll(p) = [priceEuropean(p, false, k, forward, τ, 1.0) for k in strikes]

    s = SUITE["cumulants"] = BenchmarkGroup()
    s["Heston closed form"] = @benchmarkable computeCumulants($(HestonParams(0.025, 1.0, 0.025, -0.65, 0.4)), 1.0)
    s["Heston small kappa"] = @benchmarkable computeCumulants($(HestonParams(0.025, 1e-4, 0.025, -0.65, 0.4)), 1.0)
    s["Heston small kappa, long maturity"] = @benchmarkable computeCumulants($(HestonParams(0.04, 1e-3, 0.04, -0.9, 3.0)), 30.0)

    b = SUITE["build"] = BenchmarkGroup()
    b["Cos 256"] = @benchmarkable makeCosCharFuncPricer($cf, $τ, 256, 12)
    b["Cos Junike 1e-8"] = @benchmarkable makeCosCharFuncPricer($cf, $τ, tol=1e-8)
    b["Andersen-Lake"] = @benchmarkable ALCharFuncPricer($cf)
    b["Joshi-Yang 64"] = @benchmarkable JoshiYangCharFuncPricer($cf, $τ, n=64)
    b["Adaptive Filon 1e-8"] = @benchmarkable AdaptiveFilonCharFuncPricer($cf, $τ, qTol=1e-8)
    b["FCC 1e-8"] = @benchmarkable FCCCharFuncPricer($cf, $τ, qTol=1e-8)
    b["FCC 1e-12"] = @benchmarkable FCCCharFuncPricer($cf, $τ, qTol=1e-12)

    p = SUITE["price 100 strikes"] = BenchmarkGroup()
    p["Cos 256"] = @benchmarkable $priceAll($(makeCosCharFuncPricer(cf, τ, 256, 12)))
    p["Andersen-Lake"] = @benchmarkable $priceAll($(ALCharFuncPricer(cf)))
    p["Joshi-Yang 64"] = @benchmarkable $priceAll($(JoshiYangCharFuncPricer(cf, τ, n=64)))
    p["Adaptive Filon 1e-8"] = @benchmarkable $priceAll($(AdaptiveFilonCharFuncPricer(cf, τ, qTol=1e-8)))
    p["FCC 1e-8"] = @benchmarkable priceEuropean($(FCCCharFuncPricer(cf, τ, qTol=1e-8)), $isCalls, $strikes, $forward, $τ, 1.0)
    p["FCC 1e-12"] = @benchmarkable priceEuropean($(FCCCharFuncPricer(cf, τ, qTol=1e-12)), $isCalls, $strikes, $forward, $τ, 1.0)
    p["Cos FFT 256"] = @benchmarkable priceEuropean($(CosFFTCharFuncPricer(cf, τ, 256, 12)), $isCalls, $strikes, $forward, $τ, 1.0; useNFFT=false)
end
