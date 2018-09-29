# used to make a single (discontinuous) curve that goes through a list of disconnected segments

using PyPlot
using CosineSeries
import Base.string
import Base.+

function +(a::String, b) return string(a,b) end

datapoints = [[0. 0.; 0. 0.]]
empty!(datapoints)
# paths_split has an even number of lines.
# The odd lines are the x coordinates of the segments, the even lines are the y coordinates
open("paths_split") do f
	lines = readlines(f)
	n = div(length(lines), 2)
	for i = 1:div(length(lines), 2)
		x = [parse(Float64, s) for s in split(lines[2*i-1],   ", ")]
		y = [parse(Float64, s) for s in split(lines[2*i], ", ")]
		y = 1000 .- y
		push!(datapoints, [x y])
	end
end

# Plot the segments separately.
for segment in datapoints
	x = 2 * segment[:,1]
	y = 2 * segment[:,2]
	plotCurve(x,y,10000)
end

# a single term is a_k cos(πk/N*((N-1)u + 1/2)) for u = 0..1
# Let's say we have functions f1(u), f2(u), ... for the different segments
# Assume b(t) is the standard unit box function
# For the first curve, the function is f1(u)b(u)
# For the second, it is f2(u-1)b(u-1)
# etc.
# Then the entire function is
# b(u) * ([a_11 b_11] cos(...) + [a_12 b_12] cos(...) + ...)
# + b(u-1) * ([a_21 b_21] cos(... - 1) + ...) + ...

# We can express b(x) as 1/2 (sign x + sign (1 - x)) and therefore we can
# just halve the coefficients and multiply the cosine series by (sign x + sign(1-x))
# Although we don't have to bother scaling the coefficients anyway, because we can just scale the entire
# drawing by a factor of two (only relative differences are visible)

# If we have to shift the cosine by Δ to a term of cos(kπ/N ((u - Δ) + 1/2)), we get
# cos(kπ*(N-1)/N ( u - Δ + 1 / (2N + 2)))

# We can thus see the general structure of the equation.
# Let us first define some data types.
struct CosWave # A cos(2πν(t + δ))
	A::Vector{Float64}
	ν::Float64
	δ::Float64
end

function ω(c::CosWave) return 2π*c.ν end

# Series of cosine waves that approximate a given segment
struct CosSeries
	x_min::Float64 # The value of x (t) where this function becomes active (left wall of the unit box)
	terms::Vector{CosWave}
end

# Series of disconnected
struct SegmentedSeries
	series::Vector{CosSeries}
end

function string(x::Float64)
	return string(Float32(trunc(5x))/5)
end

# "+ x" or "- x" depending on the sign of x
function strWithSgn(x)
	if x >= 0
		return "+ " + string(x)
	else
		return "- " + string(abs(x))
	end
end

function string(c :: CosWave)
	(c.ν == 0) && return string(c.A)
	return string(c.A) + " cos(" + string(ω(c)) + "t " + strWithSgn(ω(c) * c.δ) + ")"
end

# evaluates a single cosine wave in the points t
function eval(c :: CosWave, t)
	return c.A * cos.(ω(c)*(t' .+ c.δ))
end

function string(s :: CosSeries)
	repr = "(" + string(s.terms[1]) + " "
	for i = 2:length(s.terms)
		repr += " + " + string(s.terms[i])
	end
	repr += ")(sign(t " + strWithSgn(- Int64(s.x_min)) + ") - sign(t " + strWithSgn(Int64(- 1 - s.x_min)) + "))"
	return repr
end

# evaluates a cosine series in the points t
function eval(s :: CosSeries, t)
	box = (sign.(t' .- s.x_min) - sign.(t' .- s.x_min .- 1))
	bb  = [box; box]
	return sum([bb .* eval(wave, t) for wave in s.terms])
end

function string(s :: SegmentedSeries)
	repr = "[x(t), y(t)] = " + string(s.series[1])
	for i = 2:length(s.series)
		repr += " + " + string(s.series[i])
	end
	return repr + "\nt ∈ [0," + string(length(s.series)) + "]"
end

# evaluates a series of disconnected series in the points t
function eval(s :: SegmentedSeries, t)
	return sum([eval(series, t) for series in  s.series])
end

# Make the data structure for the entire equation
segmSeries = SegmentedSeries(CosSeries[])
for i = 1:length(datapoints)
	segment  = datapoints[i]
	x = segment[:,1]
	y = segment[:,2]
	N = length(x)
	a = cosAmplitudes(x)
	b = cosAmplitudes(y)
	Δ = i - 1
	series = CosSeries(Δ,[
			    CosWave([a[k+1], b[k+1]], k * (N-1)/N/2, - Δ + 1/(2N-2))
	for k=0:N-1])
	push!(segmSeries.series, series)
end

# Print the equation for the drawing
println(string(segmSeries))

# # Plot the equation whe have just found
# t  = range(0, stop=length(segmSeries.series), length=5000)
# ev = eval(segmSeries, t)
# x  = ev[1,:]
# y  = ev[2,:]
# PyPlot.plot(x, y)
