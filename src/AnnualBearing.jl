module AnnualBearing
#Mesh
G = 0.4
K = 2000
b = 200
ρ = 1
η = 0.25
S = 20
te = b*(1-η)/S/(1+ρ)
k = √(12*G/K/te^2)
params = (G=G, K=K, b=b, te=te, S=S, η=η, cs="annular")

f1(p) = msh.dcircle(p,0.,0.,b)
f2(p) = msh.dcircle(p,0.,0.,b*η)
fd(p) = msh.ddiff(f1(p),f2(p))
fh(p) = msh.huniform(p)
h0 = b/15
bbox = 1.5*[-b -b/ρ;b b/ρ]
p, t = msh.distmesh2d(fd, fh, h0, bbox)
t = msh.clean_delaunay!(p,t)
en = msh.boundary_nodes(t)

plt.figure(figsize=(3.33,3))
tplot(p,t);

#Analysis
ε = -0.1
u = meso_solver(p,t,en,(x,y)->ffun(x,y;params,ε),k)
uex = pfun.(p[:,1],p[:,2];params,ε)
error = @. abs(uex - u)/abs(uex)
error[abs.(uex) .< maximum(abs.(uex))/20] .= 0

plt.figure(figsize=(10,4))
plt.suptitle("Compression ε = $(round(ε,sigdigits=2))")
plt.subplot(1,3,1)
tplot(p,t,u)
plt.title("FEA solution")
plt.subplot(1,3,2)
tplot(p,t,uex)
plt.title("Analytical solution")
plt.subplot(1,3,3)
tplot(p,t,error)
plt.title("Error")

κy = 0.005/te
u = meso_solver(p,t,en,(x,y)->ffun(x,y;params,κy),k)
uex = pfun.(p[:,1],p[:,2];params,κy)
error = @. abs(uex - u)/abs(uex)
error[abs.(uex) .< maximum(abs.(uex))/20] .= 0

plt.figure(figsize=(10,4))
plt.suptitle("Curvature κy = $(round(κy,sigdigits=2))")
plt.subplot(1,3,1)
tplot(p,t,u,cmap="Spectral_r")
plt.title("FEA solution")
plt.subplot(1,3,2)
tplot(p,t,uex,cmap="Spectral_r")
plt.title("Analytical solution")
plt.subplot(1,3,3)
tplot(p,t,error)
plt.title("Error")

χx = 0.02/te
ω = ωfun(params)
fwx(x,y) = fwfun(x,y;params,ω)
u = meso_solver(p,t,en,(x,y)->ffun(x,y;params,χx,fwx),k)
uex = pfun.(p[:,1],p[:,2];params,χx,ω)
error = @. abs(uex - u)/abs(uex)
error[abs.(uex) .< maximum(abs.(uex))/20] .= 0

plt.figure(figsize=(10,4))
plt.suptitle("Warping rate χx = $(round(χx,sigdigits=2))")
plt.subplot(1,3,1)
tplot(p,t,u,cmap="Spectral_r")
plt.title("FEA solution")
plt.subplot(1,3,2)
tplot(p,t,uex,cmap="Spectral_r")
plt.title("Analytical solution")
plt.subplot(1,3,3)
tplot(p,t,error)
plt.title("Error");
end