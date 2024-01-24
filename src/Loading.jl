# ---------- LOADING FUNCTION ----------
#=
INPUT:
x: x coordinate
y: y coordinate
params: named tuple with keys [:G, :K, :te]
ε: axial strain
κx: curvature along x axis
κy: curvature along y axis
χx: rate of change in warping amplitude due to shear along x axis
χy: rate of change in warping amplitude due to shear along y axis
fwx: warping function due to shear along x axis
fwy: warping function due to shear along y axis
OUTPUT:
f: loading function
=#

function ffun(x::Real,y::Real;params::NamedTuple,
    ε::Real=0.,κx::Real=0.,κy::Real=0.,χx::Real=0.,χy::Real=0.,
    fwx::Function=(x,y)->1,fwy::Function=(x,y)->1)
G, K, te = params.G, params.K, params.te
f = -12*G/te^2*(ε + κx*y - κy*x - χx*fwx(x,y) + χy*fwy(x,y))
return f
end

# ---------- WARPING FUNCTION ----------
#=
INPUT:
x: x coordinate
y: y coordinate
params: named tuple with keys [:b, :ρ, :η]
ω: parameter for warping function
dir: string indicating shear direction ("x" or "y")
OUTPUT:
fw: warping function due to shear
=#

function fwfun(x::Real,y::Real;params::NamedTuple,
    ω::Real=ωfun(params),dir::AbstractString="x")
b = params.b
ρ = get(params,:ρ,1)
η = get(params,:η,0)
if dir == "y"; (x, y) = (y, x); b = b/ρ; end
if params.cs == "rectangular"
    fw = 5/6*(x^3/(2*b^2) + ω*x)
elseif params.cs == "circular" || ((params.cs == "annular") && (iszero(η)))
    r = norm([x,y])
    fw = 6/7*(r^3/(2*b^2) + ω*r)*x/r
elseif params.cs == "annular"
    r = norm([x,y])
    fw = 6/7*(r^3/(2*b^2) + ω*r - 3*η^2*b^2/(2*r))*x/r
end
return fw
end;
