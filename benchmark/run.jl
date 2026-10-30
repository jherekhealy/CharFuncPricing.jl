#julia --project=benchmark benchmark/run.jl [group]
#   group: optional name of a benchmark group, e.g. "build"
using BenchmarkTools
include(joinpath(@__DIR__, "benchmarks.jl"))

suite = isempty(ARGS) ? SUITE : SUITE[ARGS[1]]
tune!(suite)
results = run(suite, verbose=false)
for (group, res) in sort(collect(BenchmarkTools.leaves(results)), by=first)
    println(rpad(join(group, " / "), 50), BenchmarkTools.prettytime(time(median(res))), "  (", BenchmarkTools.prettymemory(memory(median(res))), ")")
end
