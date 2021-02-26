# Package CharFuncPricing

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/jherekhealy/CharFuncPricing.jl.svg?branch=master)](https://travis-ci.org/jherekhealy/CharFuncPricing.jl) | [![codecov.io](http://codecov.io/github/jherekhealy/CharFuncPricing.jl/coverage.svg?branch=master)](http://codecov.io/github/jherekhealy/CharFuncPricing.jl?branch=master) |


Julia package to provide reference European option prices for stochastic volatility models with a known characteristic function, such as the Heston stochastic volatility model. For the Heston and Schobel-Zhu models, a function provides the first, second and fourth cumulants through analytical formulas.

The code is not meant for production purpose and does not always cater for corner cases. It however supports arbitrary precision via the Nemo package.

## Installation

In a Julia REPL, enter `pkg` mode (by pressing `]`) and run:

```julia
(v1.0) pkg> add CharFuncPricing
```

[Julia](https://julialang.org) version 1.2 or higher is required.
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
### BigFloat Example

## Adaptive Flinn
### Float64 Example
### BigFloat Example

## Testing

In a Julia REPL session, enter `pkg` mode and run `test CharFuncPricing`.

Unit tests verify the call and put option prices against the [reference prices](https://financepress.com/2019/02/15/heston-model-reference-prices/) of Alan Lewis in double and arbitrary precision. In fact, the implementation here gives more precise results (minimum accuracy of 1e-60 while Alan Lewis numbers are only accurate up to 1e-25).

Cumulants are checked against a Taylor series algorithmic differentiation.

## References
Andersen, L.B.G. and Lake, M. (2018) [Robust high-precision option pricing by
Fourier transforms: Contour deformations and double-exponential quadrature](SSRN 3231626)

Fang, F. and Oosterlee, C. W. (2008) [A novel pricing method for European options based on Fourier-cosine series expansions](https://epubs.siam.org/doi/abs/10.1137/080718061)

Healy, J. (2021) [Applied Quantitative Finance for Equity Derivatives]()

Le Floc'h, F. (2018) [More Robust Pricing of European Options Based on Fourier Cosine Series Expansions](https://arxiv.org/abs/2005.13248)

Le Floc'h, F. (2020) [An adaptive Filon quadrature for stochastic volatility
models]()
