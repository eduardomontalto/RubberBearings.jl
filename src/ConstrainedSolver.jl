# ---------- CONSTRAINED SCREENED POISSON SOLVER ----------
#=
INPUT:
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
nfix: array of fixed nodes: Size: nnf x 1
f: loading function given as f(x,y) 
k: coefficient multiplying the mass matrix
ufix: array of fixed nodal values. Size: nnf x 1
ulim: limit value of the constraint such that u > ulim
OUTPUT:
u: solution vector. Size: nn x 1
=#

function meso_constrainedsolver(p::AbstractArray{Float64},t::AbstractArray{Int},nfix::Vector{Int},
    f::Function,k::Real;ufix::Vector{Float64}=zeros(length(nfix)),ulim::Real=0.)
K, M, b = meso_assembly(p,t,f)
A = K + k^2*M
Am, bm = copy(A), copy(b)
Am, bm = meso_boundaryconditions!(Am,bm,nfix;ufix)
u = Am\bm
λ = -ones(length(u))
n = 0
while true
    ncon = (1:length(u))[(u .<= ulim) .& (λ .< 0)]
    ucon = ulim*ones(length(ncon))
    Am, bm = copy(A), copy(b)
    nfc = unique([nfix;ncon])
    ufc = [ufix;ucon][1:length(nfc)]
    Am, bm = meso_boundaryconditions!(Am,bm,nfc;ufix=ufc)
    u = Am\bm
    λ = zeros(length(u))
    λ[ncon] = b[ncon] - A[ncon,:]*u 
    n += 1
    if (all(λ .<= 0)) || (n > 20); break; end
end
return u, λ
end;

