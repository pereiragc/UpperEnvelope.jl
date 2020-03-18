using UpperEnvelope
using Plots

# Example
k1 = [2.3, 2.8, 3., 4.4, 5.2]
v1 = [4.1, 1.8, 5., 4.2, 2.]

k2 = [2.3, 2.6, 3.8, 4.1, 4.8, 5.4]
# k2 = [0.5, 1.0, 1.9, 2.5]
v2 = [5.0, 0., 9.5, 3., 4.2, -2.]


pl1 = PiecewiseLinear(copy(k1), copy(v1))
pl2 = PiecewiseLinear(copy(k2), copy(v2))

# Plot example
plot(k1, v1, marker=(:black), line=(:black), label="Function 1", size=(580, 325), dpi=140)
plot!(k2, v2, marker=(:x, :darkblue), line=(:dash, :darkblue), label="Function 2")


segments = (PiecewiseLinear(copy(k1), copy(v1)), PiecewiseLinear(copy(k2), copy(v2)))
envelope=compute_envelope(segments)
plot!(envelope.xcoords, envelope.ycoords, line=(:purple, 3, :dot), label="Upper envelope")
savefig("example.png")