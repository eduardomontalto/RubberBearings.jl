module RubberBearings
using  LinearAlgebra, SparseArrays, Printf, PyCall, PyPlot
include("LinearTriangle.jl")
include("TriangularElement.jl")
include("RectangleHoleTriangleElement.jl")
include("ConstrainedSolver.jl")
include("Plotting.jl")
include("Loading.jl")
end
