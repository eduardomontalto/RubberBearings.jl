# ---------- INTEGRATE TRIANGULAR ELEMENT ----------
#= 
INPUT:
pel: 3x2 array including the (x,y) coordinates of each node
h: coordinate function given as h(x,y)
OUTPUT:
whel: element integral vector. Size: 3x1
=#
function integral_triangle(pel::AbstractArray{Float64},h::Function)
    # Number of nodes per element
    nne = 3
    # Area
    A = msh.tri_area(pel[:,1],pel[:,2])
    # Shape functions ϕᵢ(x,y) = C[1,i] + C[2,i]*x + C[3,i]*y
    V = [ones(3) pel]
    C = V\(1.0*Matrix(I,3,3))
    # Gauss integration points
    xg = sum(pel[:,1])/6 .+ pel[:,1]/2
    yg = sum(pel[:,2])/6 .+ pel[:,2]/2
    wg = A/3*ones(3)
    whel = zeros(nne,1)
    for i in 1:nne, g in eachindex(xg)
        whel[i] += (C[:,i]'*[1,xg[g],yg[g]])*h(xg[g],yg[g])*wg[g]
    end
    return whel
end

# ---------- ASSEMBLE INTEGRATION VECTOR ---------
#=
INPUT:
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
h: coordinate function given as h(x,y) 
OUTPUT:
wh: mesh integral vector. Size: nn x 1
=#
function integral_assembly(p::AbstractArray{Float64},t::AbstractArray{Int},h::Function)
    nn = size(p,1)
    nel = size(t,1)
    wh = zeros(nn)
    for i in 1:nel
        pel = [p[t[i,:],1] p[t[i,:],2]]
        wh[t[i,:]] += integral_triangle(pel,h)
    end
    return wh
end

# ---------- INTEGRATE FEA SOLUTION ----------
#=
u: array of FEA nodal values. Size nn x 1
h: coordinate function given as h(x,y)
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
OUTPUT:
int: integral of u(x,y)*h(x,y)
=#
function integrate(u::Vector{Float64},h::Function,p::AbstractArray{Float64},t::AbstractArray{Int})
    int = u'*integral_assembly(p,t,h)
    return int
end;

# ---------- GENERALIZED STRESS CALCULATION ----------
#=
u: array of pressure nodal values. Size nn x 1
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
OUTPUT:
N, Mx, My: generalized stress
=#
function Nfun(u::Vector{Float64},p::AbstractArray{Float64},t::AbstractArray{Int})
    N = integrate(-u,(x,y)->1,p,t)
    return N
end;

function Mxfun(u::Vector{Float64},p::AbstractArray{Float64},t::AbstractArray{Int})
    Mx = integrate(-u,(x,y)->y,p,t)
    return Mx
end;

function Myfun(u::Vector{Float64},p::AbstractArray{Float64},t::AbstractArray{Int})
    My = integrate(-u,(x,y)->-x,p,t)
    return My
end

function Afun(λ::Vector{Float64},p::AbstractArray{Float64},t::AbstractArray{Int})
    a = ones(length(λ))
    a[λ .!= 0.] .*= 0.
    A = integrate(a,(x,y)->1,p,t)
    return A
end;
#Mesoscale Model
# ---------- MESOSCALE MODEL FOR 2D DEFORMATION ----------
#=
e: vector of generalized strains [ε, γx, κy]
params: named tuple with keys [:G, :K, :te]
Mesh: named tuple with keys [:p, :t, :en] 
    p: array of (x,y) coordinates of cross section. Size: nn x 2
    t: array of connectivity of cross section. Size: nel x 3
    en: array of nodes in the boundary of cross section. Size: nnb x 1
Δe: vector of perturbation of generalized strains [Δε, Δγx, Δκy] 
OUTPUT:
s: vector of generalized stresses [N, My, Vx]. Size: 3 x 1
ks: stiffness matrix corresponding to ∂s/∂e. Size: 3 x 3
=#
function meso_2d_nw(e::Vector{Float64},params::NamedTuple,Mesh::NamedTuple;
    Δe::Vector{Float64}=[0.0001,0.0001,0.00001])
G, K, te = params.G, params.K, params.te
ε, γx, κy = e
Δε, Δγx, Δκy = Δe
k = √(12*G/K/te^2)
u, λ = meso_constrainedsolver(Mesh.p,Mesh.t,Mesh.en,(x,y)->ffun(x,y;params,ε,κy),k)
s = [Nfun(u,Mesh.p,Mesh.t),G*Afun(λ,Mesh.p,Mesh.t)*γx,Myfun(u,Mesh.p,Mesh.t)]
ks = zeros(3,3)
ks[2,2] += G*Afun(λ,Mesh.p,Mesh.t)
u, λ = meso_constrainedsolver(Mesh.p,Mesh.t,Mesh.en,(x,y)->ffun(x,y;params,ε=ε-Δε,κy),k)
s1 = [Nfun(u,Mesh.p,Mesh.t),G*Afun(λ,Mesh.p,Mesh.t)*γx,Myfun(u,Mesh.p,Mesh.t)]
u, λ = meso_constrainedsolver(Mesh.p,Mesh.t,Mesh.en,(x,y)->ffun(x,y;params,ε=ε+Δε,κy),k)
s2 = [Nfun(u,Mesh.p,Mesh.t),G*Afun(λ,Mesh.p,Mesh.t)*γx,Myfun(u,Mesh.p,Mesh.t)]
u, λ = meso_constrainedsolver(Mesh.p,Mesh.t,Mesh.en,(x,y)->ffun(x,y;params,ε,κy=κy-Δκy),k)
s3 = [Nfun(u,Mesh.p,Mesh.t),G*Afun(λ,Mesh.p,Mesh.t)*γx,Myfun(u,Mesh.p,Mesh.t)]
u, λ = meso_constrainedsolver(Mesh.p,Mesh.t,Mesh.en,(x,y)->ffun(x,y;params,ε,κy=κy+Δκy),k)
s4 = [Nfun(u,Mesh.p,Mesh.t),G*Afun(λ,Mesh.p,Mesh.t)*γx,Myfun(u,Mesh.p,Mesh.t)]
ks[:,1] = (s2 - s1)./(2*Δε)
ks[:,3] = (s4 - s3)./(2*Δκy)
return s, ks
end;
