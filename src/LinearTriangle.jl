module LinearTriangle

# Write your package code here.
# ----------LINEAR TRIANGULAR ELEMENT----------
#= 
INPUT:
pel: 3x2 array including the (x,y) coordinates of each node
f: loading function given as f(x,y)
OUTPUT:
Kel: element stiffness matrix. Size: 3x3
Mel: element mass matrix. Size: 3x3
bel: element loading vector. Size: 3x1
=#
function linear_triangle(pel::AbstractArray{Float64},f::Function)
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
    # Stiffness matrix
    Kel = zeros(nne,nne)
    for i in 1:nne, j in i:nne
        Kel[i,j] = C[2:3,i]'*C[2:3,j]*A
    end
    Kel = Symmetric(Kel)
    # Mass matrix
    Mel = zeros(nne,nne)
    for i in 1:nne, j in i:nne, g in eachindex(xg)
        Mel[i,j] += (C[:,i]'*[1,xg[g],yg[g]])*(C[:,j]'*[1,xg[g],yg[g]])*wg[g]
    end
    Mel = Symmetric(Mel)
    # Loading vector
    bel = zeros(nne,1)
    for i in 1:nne, g in eachindex(xg)
        bel[i] += (C[:,i]'*[1,xg[g],yg[g]])*f(xg[g],yg[g])*wg[g]
    end
    return Kel, Mel, bel
end

# ----------ASSEMBLE MATRICES AND LOAD VECTOR---------
#=
INPUT:
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
f: loading function given as f(x,y) 
OUTPUT:
K: stiffness matrix. Size: nn x nn
M: mass matrix. Size: nn x nn
b: loading vector. Size: nn x 1
=#
function meso_assembly(p::AbstractArray{Float64},t::AbstractArray{Int},f::Function)
    nn = size(p,1)
    nel = size(t,1)
    nne = size(t,2)
    K = Tuple{Int,Int,Real}[]
    M = Tuple{Int,Int,Real}[]
    b = zeros(nn)
    for i in 1:nel
        pel = [p[t[i,:],1] p[t[i,:],2]]
        Kel, Mel, bel = linear_triangle(pel,f)
        for j in 1:nne, k in 1:nne
            push!(K, (t[i,j], t[i,k], Kel[j,k]))
            push!(M, (t[i,j], t[i,k], Mel[j,k]))
        end
        b[t[i,:]] += bel
    end
    K = sparse((x->x[1]).(K), (x->x[2]).(K), (x->x[3]).(K))
    M = sparse((x->x[1]).(M), (x->x[2]).(M), (x->x[3]).(M))
    return K, M, b
end

# ---------- APPLY BOUNDARY CONDITIONS ---------
#=
INPUT:
A: matrix in linear systems of eqs. (K + k^2 M). Size: nn x nn
b: loading vector. Size: nn x 1
nfix: array of fixed nodes: Size: nnf x 1
ufix: array of fixed nodal values. Size: nnf x 1
OUTPUT:
A: matrix in linear systems of eqs. (K + k^2 M). Size: nn x nn
b: loading vector. Size: nn x 1
=#
function meso_boundaryconditions!(A::AbstractArray{Float64},b::Vector{Float64},nfix::Vector{Int};
    ufix::Vector{Float64}=zeros(length(nfix)))
    b -= A[:,nfix]*ufix
    A[:,nfix] .= 0.
    A[nfix,:] .= 0.
    A[diagind(A)[nfix]] .= 1.
    b[nfix] = ufix
    dropzeros!(A)
    return A, b
end;

# ---------- SCREENED POISSON SOLVER ----------
#=
INPUT:
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
nfix: array of fixed nodes: Size: nnf x 1
f: loading function given as f(x,y) 
k: coefficient multiplying the mass matrix
ufix: array of fixed nodal values. Size: nnf x 1
OUTPUT:
u: solution vector. Size: nn x 1
=#

function meso_solver(p::AbstractArray{Float64},t::AbstractArray{Int},nfix::Vector{Int},
        f::Function,k::Real;ufix::Vector{Float64}=zeros(length(nfix)))
    K, M, b = meso_assembly(p,t,f)
    A = K + k^2*M
    A, b = meso_boundaryconditions!(A,b,nfix;ufix)
    u = A\b
    return u
end;
end
