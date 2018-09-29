# Finds an interpolating cosine series parametrisation for a given set of x and y values
# The numer of terms in the cosine series is exactly equal to the number of data points
module CosineSeries

using PyPlot
using FFTW
export cosAmplitudes, cosSeries, plotCurve

# Formulas for the forward and backward DCT:
# X_k = \sqrt{2/N} sum_{n=0}^{N-1} x_n cos(π/N (n+1/2) k)
# x_n = \sqrt{2/N} [ 1/sqrt(2) X_0 + sum_{k=0}^{N-1} X_k cos(π/2 (n+1/2) k) ]
# assuming indices start from 0

# This returns the coefficients to an interpolating cosine series for a sequence of sampled points
# If the result is the vector a, then the following function will interpolate x:
# f(t) = sum_{k=0}^{N-1} a_k cos(πk/N (t+1/2)) with t = 0..N-1
function cosAmplitudes(x)
	N = length(x)
	# Apply the formula for the inverse DCT and substitue n=t (discrete to continuous variable)
	amps = dct(x)
	amps[1] /= sqrt(2)
	amps    *= sqrt(2/N)
	return amps
end

# A function f(t) that interpolates x
function cosSeries(x)
	N = length(x)
	amps = cosAmplitudes(x)
	f(t) = sum([amps[k+1] * cos.(π*k/N * (t.+1/2)) for k=0:N-1])
	return f
end

# Plots the interpolating cosine curve given the x and y values of the sampled points
function plotCurve(x, y, resolution)
	N = length(x)
	N != length(y) && throw(ArgumentError("x and y must be the same size"))
	t = range(0, stop=N-1, length=resolution)
	x_curve = cosSeries(x)
	y_curve = cosSeries(y)
	PyPlot.plot(x_curve(t), y_curve(t), "b")
end

end # module
