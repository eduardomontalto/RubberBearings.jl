module RectangularAndHole
G = 0.4
K = 2000
b = 200
ρ = 1
η = 0.2
l = b/ρ
a = b*η
S = 20
te = (4*b*l - π*a^2)/(4*(b+l) + 2*π*a)/S
k = √(12*G/K/te^2)
params = (G=G, K=K, te=te)

f1(p) = msh.drectangle(p,-b,b,-b/ρ,b/ρ)
f2(p) = msh.dcircle(p,0.,0.,b*η)
fd(p) = msh.ddiff(f1(p),f2(p))
fh(p) = msh.huniform(p)
h0 = b/15
bbox = 1.5*[-b -b/ρ;b b/ρ]
pfix = [-b -b/ρ;b -b/ρ;b b/ρ;-b b/ρ]
p, t = msh.distmesh2d(fd, fh, h0, bbox, pfix=pfix)
t = msh.clean_delaunay!(p,t)
en = msh.boundary_nodes(t)

plt.figure(figsize=(3.33,3))
tplot(p,t);
#Analysis
params = merge(params,(S=S,))
β = βfun(params)/(1+ρ)
ε = -0.00513*3
κx = 0.
κy = -ε/b/((0.4-0.05*β)*ρ+1)/((45/4/β^2 + 45/28)^(-1) + 1/3)*0.77
u, λ = meso_constrainedsolver(p,t,en,(x,y)->ffun(x,y;params,ε,κx,κy),k)

plt.figure(figsize=(6.67,4))
plt.suptitle("Strains: ε = $(round(ε,sigdigits=2)), κx = $(round(κx,sigdigits=2)), κy = $(round(κy,sigdigits=2))")
plt.subplot(1,2,1)
tplot(p,t,u)
plt.title("Pressure")
plt.subplot(1,2,2)
tplot(p,t,λ,cmap="YlGnBu_r")
plt.title("Lagrange multipliers")
# @printf("N = %.2e kN \n",Nfun(u,p,t)/1e3)
# @printf("Mx = %.2e kN-m \n",Mxfun(u,p,t)/1e6)
# @printf("My = %.2e kN-m \n",Myfun(u,p,t)/1e6)
# @printf("Contact = %.2f \n",Afun(λ,p,t)/(4*b*l - π*a^2));
end