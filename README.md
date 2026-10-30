# Package CharFuncPricing

| Status | Coverage |
| :----: | :----: |
| ![Build Status](https://github.com/jherekhealy/CharFuncPricing.jl/actions/workflows/julia-runtests.yml/badge.svg) | [![codecov.io](https://codecov.io/github/jherekhealy/CharFuncPricing.jl/coverage.svg?branch=master)](https://codecov.io/github/jherekhealy/CharFuncPricing.jl?branch=master) |

Julia package to provide reference European option prices for stochastic volatility models with a known characteristic function, such as the Heston stochastic volatility model.

The code is not meant for production purpose and does not always cater for corner cases. It however supports arbitrary precision via BigFloat or the Nemo package, and automatic differentiation via ForwardDiff.

### Models
| Model | Parameters | Characteristic function |
| :---- | :---- | :---- |
| Heston | `HestonParams(v0, κ, θ, ρ, σ)` | `DefaultCharFunc(params)` |
| Heston with piecewise-constant parameters | `HestonTSParams` | `DefaultCharFunc(params)` |
| Double Heston | `DoubleHestonParams(heston1, heston2)` | `DefaultCharFunc(params)` |
| Schöbel-Zhu (and with term-structure) | `SchobelZhuParams`, `SchobelZhuTSParams` | `DefaultCharFunc(params)` |
| Stochastic volatility with correlated jumps | `SVCJParams` | `DefaultCharFunc(params)` |
| CGMY, Variance Gamma | `CGMYParams`, `VarianceGammaParams` | `DefaultCharFunc(params)` |
| Rough Heston | `RoughHestonParams(H, ρ, ν, varianceCurve)` | `PadeCharFunc(params)`, `AdamsCharFunc(params, n)` |

For the Heston and Schöbel-Zhu models, `computeCumulants(params, τ)` provides the first, second and fourth cumulants through analytical formulas (for Heston with small `κτ`, a numerically stable ODE solution is used instead of the closed form).

### Pricers
| Method | Constructor |
| :---- | :---- |
| Cos (Fang-Oosterlee), cumulants truncation | `makeCosCharFuncPricer(cf, τ, m, l)` |
| Cos on an explicit interval `[a, b]` | `makeCosCharFuncPricer(cf, τ, m, a, b)` |
| Cos with automatic truncation and number of terms (Junike) | `makeCosCharFuncPricer(cf, τ, tol=1e-8)` |
| Andersen-Lake contour integration | `ALCharFuncPricer(cf)` |
| Adaptive Flinn / Filon quadratures | `AdaptiveFlinnCharFuncPricer(cf, τ)`, `AdaptiveFilonCharFuncPricer(cf, τ)` |
| Joshi-Yang with Black-Scholes control variate | `JoshiYangCharFuncPricer(cf, τ)` |
| Adaptive Filon-Clenshaw-Curtis, efficient for many strikes | `FCCCharFuncPricer(cf, τ, qTol=1e-10)` |
| Swift (Shannon wavelets) | `makeSwiftCharFuncPricer(cf, τ, m, k)` |
| Cos variants: FFT/NFFT on many strikes, Lipton, damped (α), Gil-Pelaez | `CosFFTCharFuncPricer`, `makeCosLiptonCharFuncPricer`, `makeCosAlphaCharFuncPricer`, `makeCosGPCharFuncPricer` |

All pricers are used through `priceEuropean(pricer, isCall, strike, forward, τ, discountDf)`, or `priceEuropean(pricer, isCalls, strikes, forward, τ, discountDf)` for many strikes at once.

## Installation

In a Julia REPL, enter `pkg` mode (by pressing `]`) and run:

```julia
pkg> add CharFuncPricing
```

[Julia](https://julialang.org) version 1.10 or higher is required.

## Cos method
### Float64 Usage

Start by creating a `HestonParams` structure, which represents the parameters of the Heston model: v0, κ, θ, ρ, σ.

```julia
params = HestonParams(v0, κ, θ, ρ, σ)
```

Then make a `CosCharFuncPricer` structure via `makeCosCharFuncPricer`. This function will store the relevant `m` values of the characteristic function for the range [a,b] defined by `l` deviations using the cumulants rule `a = c1 - l * sqrt(c2+sqrt(c4))`, `b = c1 + l * sqrt(c2+sqrt(c4))`.
```julia
cf = DefaultCharFunc(params)
pricer = makeCosCharFuncPricer(cf, τ, m, l)
```

And price vanilla call and puts of a given strike as following
```julia
priceEuropean(pricer, false, strike, forward, τ, df)
```
The second parameter specifies whether we want to price a call (true) or a put (false). The last parameter specifies the discount factor to maturity.

The first, second and fourth cumulants are given by
```julia
c1,c2,c4 = computeCumulants(params, τ)
```

Alternatively, the truncation range and the number of terms may be determined automatically for a given tolerance, following Junike (2024):
```julia
pricer = makeCosCharFuncPricer(cf, τ, tol=1e-8)
```

### Nemo Usage
The only difference is to make sure the parameters are of ArbField type (real arbitrary precision), the function `makeCosCharFuncPricer` should also be called on the AcbField type.
```julia
using Nemo
R = ArbField(256)
CC = AcbField(256)
cf = NemoCharFunc{HestonParams{arb}}(params, CC)
pricer = makeCosCharFuncPricer(cf, τ, m, l)
```

### Float64 Example
Here is how to price a put option with parameters
```julia
r=0.01; q=0.02
κ=4.0; θ=0.25; σ=1.0; ρ=-0.5; v0=0.04
τ = 1.0
spot = 100.0; strike = 80.0
spot *= exp((r - q) * τ)
df = exp(-r * τ)
params = HestonParams(v0, κ, θ, ρ, σ)
cf = DefaultCharFunc(params)
l = 32; m = 1024
pricer = makeCosCharFuncPricer(cf, τ, m, l)
priceEuropean(pricer, false, strike, spot, τ, df)
```

The result is `7.95887811325676`.

### Nemo Example
```julia
using Nemo
using CharFuncPricing

R = ArbField(256)
CC = AcbField(256)
r = R(QQ(1, 100)); q = R(QQ(2, 100))
κ = R(4.0); θ = R(QQ(1, 4)); σ = R(1.0); v0 = R(QQ(4, 100)); ρ = R(-0.5);
τ = R(1.0)
spot = R(100); strike = R(80)
isCall = false
spot *= exp((r - q) * τ)
df = exp(-r * τ)
params = HestonParams(v0, κ, θ, ρ, σ)
cf = NemoCharFunc{HestonParams{arb}}(params, CC)
l = 48; m = 1024*4
pricer = makeCosCharFuncPricer(cf, τ, m, l)
priceEuropean(pricer, isCall, strike, spot, τ, df)
```

The result is:

`7.95887811325676828521326060761429303089865693725960319205094095681790030 +/- 4.83e-72`.

With `l=64; m=1024*8`, we obtain:

`7.95887811325676828521326060761429303089865693725960319205094095681878397 +/- 3.71e-72`.


## Andersen-Lake
### Float64 Example
It starts similarly as for the Cos method, the only change is how to build the pricer variable.
```julia
pricer = ALCharFuncPricer(cf)
priceEuropean(pricer, false, strike, spot, τ, df)
```

The result is `7.9588781132567705`.

### BigFloat Example
We use BigFloat for the Heston parameters and option characteristics, as well as in the quadrature tolerance.
```julia
r=BigFloat("0.01"); q=BigFloat("0.02")
κ=BigFloat(4.0); θ=BigFloat("0.25"); σ=BigFloat(1.0); ρ=BigFloat("-0.5"); v0=BigFloat("0.04")
τ = BigFloat(1.0)
spot = BigFloat(100.0); strike = BigFloat(80.0)
spot *= exp((r - q) * τ)
df = exp(-r * τ)
params = HestonParams(v0, κ, θ, ρ, σ)
cf = DefaultCharFunc{HestonParams{BigFloat},Complex{BigFloat}}(params)
quad = TanhSinhQuadrature(800, BigFloat(1e-200))
pricer = ALCharFuncPricer(cf,quad)
priceEuropean(pricer, false, strike, spot, τ, df)
```

The result is

`7.95887811325676828521326060761429303089865693725960319205094095681918541918632`

## Adaptive Flinn
This is the adaptive Flinn quadrature using the transformation to (-1, 1) interval. No truncation is involved.

### Float64 Example
With a quadrature tolerance of 1e-8:
```julia
pricer = AdaptiveFlinnCharFuncPricer(cf, τ, qTol = 1e-8)
priceEuropean(pricer, false, strike, spot, τ, df)
```
The result is `7.958878080086734`.

### BigFloat Example
The adaptive Flinn pricer works with high accuracy, but does not perform very well then. It is more intended for the calculation of prices with a absolute error tolerance of around 1e-8 or 1e-10.
```julia
pricer = AdaptiveFlinnCharFuncPricer(cf, τ, qTol = BigFloat(1e-24))
priceEuropean(pricer, false, strike, spot, τ, df)
```
The result is

`7.958878113256768285213257572750089190600415520655780637746847607110529890012863`

  and the effective accuracy is 3e-24. In total, `length(pricer.kcos[1,:])=31511` points are used, compared to 295 for an tolerance of 1e-10. For a tolerance of 1e-32, 450799 points are used. This means that algorithm is asymptotically linear on this example.


## Joshi-Yang integration
A direct implementation of Joshi-Yang pricing with a Black-Scholes control variate integrated along the real axis combined with a Gauss-Laguerre integration is given.

```julia
using CharFuncPricing
cf = DefaultCharFunc(HestonParams(0.04, 0.5, 0.04, -0.9, 1.0))
T= 10.0
pricer = JoshiYangCharFuncPricer(cf,T, n=10)
strikes = [8000.0, 10000.0, 12000.0]
prices = map(x->priceEuropean(pricer,x >= 10000.0, x, 10000.0,T,1.0) ,strikes)
```
The output reads
`772.0090642032833
1308.1993158318369
290.12673794349496
`
which is in-line with the prices reported in Table 4.2 of the paper.

While it is accurate at or around the money in general, for a small number of integration points, this is not always accurate for long-term deals. Furthermore, for long-term deals, the control variate volatility may *"explode"*. A counter-example is  `T=30.0, HestonParams{Float64}(1.0, 0.5, 1.0, 0.95, 1.0)`.

## Filon-Clenshaw-Curtis
`FCCCharFuncPricer` integrates the Lewis formula with an adaptive Filon-Clenshaw-Curtis quadrature. The integrand amplitude `φ(u - i/2)/(u² + 1/4)` is approximated on `[0, b]` by piecewise Chebyshev polynomials of degree `N` (16 by default), on panels bisected until the last Chebyshev coefficients fall below the tolerance. The panels do not depend on the strike: the oscillatory factor is integrated exactly against the polynomials (through a stable recurrence for the moments, or Gauss-Legendre where the oscillation is weak). The characteristic function is thus evaluated only when the pricer is built, and pricing an additional strike is cheap.

```julia
pricer = FCCCharFuncPricer(cf, τ, qTol = 1e-10)
priceEuropean(pricer, false, strike, spot, τ, df)
```
On the Heston example of Alan Lewis, the error is below 1e-13 for a tolerance of 1e-8, while the adaptive Filon quadrature with the same tolerance has an error of 2.5e-6. The truncation `b` is the power of two beyond which the integrand is negligible, it may be given explicitly with the keyword `b`.

As for the other pricers based on the Lewis formula (`α = -1/2`), it fails when `φ(-i/2)` is not finite, for example with the Variance Gamma parameters `(2.0, 1.0, 0.5)`, where the Cos method still works.

## Choosing a pricer by name
`makeCharFuncPricer(method, cf, τ)` builds a pricer with settings suited to fast calibration, where `method` is one of `:AndersenLake`, `:JoshiYang`, `:Flinn`, `:FlinnTransformed`, `:Swift`, `:Cos`, `:Cos128`, `:CosJunike`, `:FCC` (a string such as `"Joshi-Yang"` is also accepted). Other packages may add their own methods by extending `makeCharFuncPricer(::Val{:Name}, cf, τ)`.

```julia
pricer = makeCharFuncPricer(:JoshiYang, cf, τ)
priceEuropean(pricer, false, strike, spot, τ, df)
```
The settings favour speed over accuracy: on the Float64 example above, `:JoshiYang` and `:Cos` are accurate to 1e-13, while `:Flinn` (with a time tolerance of 1e-4) is only accurate to 5e-4.

## Rough Heston
The rough Heston model of El Euch and Rosenbaum is parameterized by the Hurst exponent `H`, the correlation `ρ`, the vol-of-vol `ν` and the forward variance curve, a function of time. Two characteristic functions are available:
- `PadeCharFunc(params)`, the Padé (3,3) approximation of Gatheral and Radoičić, integrated in time with a Fejér quadrature (`PadeCharFunc{Complex}(params, FejerQuadrature(64))` for more points). It requires `H < 1/2`.
- `AdamsCharFunc(params, n)`, the fractional Adams scheme of Diethelm, Ford and Freed on `n` time steps, slower but converging to the exact solution, with an error of order `1 + α = H + 3/2` in the time step. `AdamsCharFunc(params, n, richardson=true)` extrapolates from `n` and `2n` steps. The explicit scheme is not finite for large `u` (beyond `u ≈ 7.5 √n` for the example below): the FCC and Joshi-Yang pricers handle this, the Cos method and the adaptive Filon quadrature do not.

The cumulants are not available, the Cos method thus requires an explicit truncation range:
```julia
using CharFuncPricing
params = RoughHestonParams(0.1, -0.65, 0.3, t -> 0.025)
pricer = makeCosCharFuncPricer(PadeCharFunc(params), 1.0, 256, -6.0, 6.0)
prices = [priceEuropean(pricer, true, k, 1.0, 1.0, 1.0) for k in [0.9, 1.0, 1.1]]
```
The result is `[0.12391243781049857, 0.04985891958466849, 0.010398774983298054]`, while the Adams scheme with 500 steps leads to `[0.1237623825975548, 0.049830724715183274, 0.010657003023154768]`. The characteristic function is expensive to evaluate, and `FCCCharFuncPricer(AdamsCharFunc(params, 500), 1.0, qTol=1e-10)` needs fewer evaluations for an accurate price, `[0.1237596731157109, 0.0498327919993613, 0.010664364715242014]`, without a truncation range to choose.

`RoughHestonCVCharFunc(cf)` adds a Black-Scholes control variate, to be used with the Fourier integral pricers such as `AdaptiveFilonCharFuncPricer`.

## Stochastic volatility with correlated jumps (SVCJ)
A straightforward (but not necessarily efficient) implementation of the SVCJ characteristic function of Duffie, Pan & Singleton is available and may be used with Cos, JoshiYang or AdaptiveFilon pricers.

```julia
using CharFuncPricing
hestonParams = HestonParams(0.007569, 3.46, 0.008, -0.82, 0.14)
params = SVCJParams(hestonParams, 0.47, -0.1, 0.0001^2, 0.05, -0.38) #from Broadie-Kaya paper
df = exp(-0.0319)
cf = DefaultCharFunc(params)
pricerCos = CharFuncPricing.makeCosCharFuncPricer(cf,1.0,1024,16)
price = priceEuropean(pricerCos,true, 100.0, 100.0/df,1.0,df)
```

The result reads `6.861875621424775` while the reference from Broadie & Kaya is 6.8619.

## Benchmarks
The folder `benchmark` contains a [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl) suite (`SUITE` in `benchmark/benchmarks.jl`, compatible with PkgBenchmark) for the cumulants, the construction of the pricers and the pricing of 100 strikes. Run it with
```
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/run.jl          # all groups
julia --project=benchmark benchmark/run.jl build    # a single group
```

## Testing

In a Julia REPL session, enter `pkg` mode and run `test CharFuncPricing`.

Unit tests verify the call and put option prices against the [reference prices](https://financepress.com/2019/02/15/heston-model-reference-prices/) of Alan Lewis in double and arbitrary precision. In fact, the implementation here gives more precise results (minimum accuracy of 1e-60 while Alan Lewis numbers are only accurate up to 1e-25).

Cumulants are checked against a Taylor series algorithmic differentiation. The rough Heston characteristic functions are checked against the Heston limit `H = 1/2` and against each other.

## References
Andersen, L.B.G. and Lake, M. (2018) "Robust high-precision option pricing by
Fourier transforms: Contour deformations and double-exponential quadrature"(SSRN 3231626)

Broadie, M. and Kaya, O. (2006) [Exact simulation of stochastic volatility and other affine jump diffusion processes](https://www.researchgate.net/publication/220243902_Exact_Simulation_of_Stochastic_Volatility_and_Other_Affine_Jump_Diffusion_Processes)

Diethelm, K., Ford, N. J. and Freed, A. D. (2002) "A predictor-corrector approach for the numerical solution of fractional differential equations", Nonlinear Dynamics 29

Duffie, D. Pan, J. and Singleton, K.(2000) [Transform Analysis and Asset Pricing for Affine Jump-Diffusion](https://www.darrellduffie.com/uploads/pubs/DuffiePanSingleton2000.pdf)

El Euch, O. and Rosenbaum, M. (2019) "The characteristic function of rough Heston models", Mathematical Finance 29(1)

Fang, F. and Oosterlee, C. W. (2008) [A novel pricing method for European options based on Fourier-cosine series expansions](https://epubs.siam.org/doi/abs/10.1137/080718061)

Gatheral, J. and Radoičić, R. (2019) "Rational approximation of the rough Heston solution", International Journal of Theoretical and Applied Finance 22(3)

Healy, J. (2021) [Applied Quantitative Finance for Equity Derivatives](https://jherekhealy.github.io/)

Joshi, M. and Yang, C. (2011) [Fourier transforms, option pricing and controls](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1941464)

Junike, G. (2024) "On the number of terms in the COS method for European option pricing", Numerische Mathematik

Le Floc'h, F. (2018) [More Robust Pricing of European Options Based on Fourier Cosine Series Expansions](https://arxiv.org/abs/2005.13248)

Le Floc'h, F. (2020) "An adaptive Filon quadrature for stochastic volatility models"

Domínguez, V., Graham, I. G. and Smyshlyaev, V. P. (2011) "Stability and error estimates for Filon-Clenshaw-Curtis rules for highly oscillatory integrals", IMA Journal of Numerical Analysis 31(4)

Ortiz-Gracia, L. and Oosterlee, C. W. (2016) "A highly efficient Shannon wavelet inverse Fourier technique for pricing European options", SIAM Journal on Scientific Computing 38(1)
