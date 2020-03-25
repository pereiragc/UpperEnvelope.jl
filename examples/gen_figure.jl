using BenchmarkTools
using UpperEnvelope
using Plots
import Random

Random.seed!(42);

tmp1 = sortslices(reshape(rand(15), (5, 3)), dims=1)
k1=tmp1[:, 1]
v1=tmp1[:, 2]
w1=1 .+ tmp1[:, 3]

tmp2 = sortslices(reshape(rand(15), (5, 3)), dims=1)
k2=tmp2[:, 1]
v2=tmp2[:, 2]
w2=1 .+ tmp2[:, 3]

envelope=compute_envelope((k1, v1, w1), (k2, v2, w2), true)

# Plot example (too lazy to make code prettier)
p0=plot(k1, [v1 w1], marker=(:black), line=(:black), label="Function 1")
p1=plot(k1, [v1 w1], marker=(:black), line=(:black), label="Function 1")
plot!(p0, k2, [v2 w2], marker=(:x, :darkblue), line=(:dash, :darkblue), label="Function 2")
plot!(p1, k2, [v2 w2], marker=(:x, :darkblue), line=(:dash, :darkblue), label="Function 2")
plot!(p1,envelope[1], [envelope[2] envelope[3]], line=(:purple, 3), marker=(:diamond, :purple, 3), label="Upper envelope")
plot(p0, p1, size=0.8.*(980, 440), dpi=190)

# savefig("example.png")


