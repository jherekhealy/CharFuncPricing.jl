# Package CharFuncPricing

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/jherekhealy/CharFuncPricing.jl.svg?branch=master)](https://travis-ci.org/jherekhealy/CharFuncPricing.jl) | [![codecov.io](http://codecov.io/github/jherekhealy/CharFuncPricing.jl/coverage.svg?branch=master)](http://codecov.io/github/jherekhealy/CharFuncPricing.jl?branch=master) |


Julia package to provide reference European option prices for stochastic volatility models with a known characteristic function, such as the Heston stochastic volatility model.

 The code is not meant for production purpose and does not cater for corner cases. It however supports arbitrary precision via the Nemo package.

## Installation

In a Julia REPL, enter `pkg` mode (by pressing `]`) and run:

```julia
(v1.0) pkg> add CharFuncPricing
```

[Julia](https://julialang.org) version 1.2 or higher is required.

## Float64 Usage

Start by creating a `HestonParams` structure, which represents the parameters of the Heston model: v0, κ, θ, ρ, σ.

```julia
params = HestonParams(v0, κ, θ, ρ, σ)
```

Then make a `CosCharFuncPricer` structure via `makeCosCharFuncPricer`. This function will store the relevant `m` values of the characteristic function for the range [a,b] defined by `l` deviations using the cumulants rule a = c_1 - l*sqrt(c_2+sqrt(c_4)), b = c_1 + l*sqrt(c_2+sqrt(c_4)).
```julia
pricer = makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)  
```

And price vanilla call and puts of a given strike as following
```julia
priceEuropean(pricer, false, strike, forward, df)
```
The second parameter specifies whether we want to price a call (true) or a put (false). The last parameter specifies the discount factor to maturity.

## Nemo Usage
## Examples


```julia
using CharFuncPricing

R = ArbField(256)
CC = AcbField(256)
r = R(QQ(1, 100))
q = R(QQ(2, 100))
κ = R(4.0)
θ = R(QQ(1, 4))
σ = R(1.0)
v0 = R(QQ(4, 100))
ρ = R(-0.5)
τ = R(1.0)
spot = R(100)
isCall = false
spot *= exp((r - q) * τ)
df = exp(-r * τ)
params = HestonParams(v0, κ, θ, ρ, σ)
l = 48
m = 1024*4
pricer = makeCosCharFuncPricer(CC, R, const_pi(R), params, τ, m, l)
priceEuropean(pricer, isCall, strike, spot, df)
```

### Plot

## Testing

In a Julia REPL session, enter `pkg` mode and run `test CharFuncPricing`.
