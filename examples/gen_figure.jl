using UpperEnvelope
using Plots

# Example
k1 = [2.3, 2.8, 3., 4.4, 5.2]
v1 = [4.1, 0.5, 5., 4.2, 2.]

k2 = [2.3, 2.6, 3.8, 4.1, 4.8, 5.4]
v2 = [5.0, 0., 9.5, 3., 4.2, -2.]

segments = (PiecewiseLinear(copy(k1), copy(v1)), PiecewiseLinear(copy(k2), copy(v2)))
envelope = compute_envelope(segments)

# Plot example (too lazy to make code prettier)
p0=plot(k1, v1, marker=(:black), line=(:black), label="Function 1")
p1=plot(k1, v1, marker=(:black), line=(:black), label="Function 1")
plot!(p0, k2, v2, marker=(:x, :darkblue), line=(:dash, :darkblue), label="Function 2")
plot!(p1, k2, v2, marker=(:x, :darkblue), line=(:dash, :darkblue), label="Function 2")
plot!(p1,envelope.xcoords, envelope.ycoords, line=(:purple, 3), marker=(:diamond, :purple, 3), label="Upper envelope")
plot(p0, p1, size=0.8.*(980, 440), dpi=190)

savefig("example.png")
