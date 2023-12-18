module RectangeHoleTriangleElement
#Mesh
G = 0.4
K = 2000
b = 200
ρ = 1
η = 0.2
l = b/ρ
a = b*η
S = 20
A = (4*b*l - π*a^2)
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
Mesh = (p=p, t=t, en=en)

plt.figure(figsize=(3.33,3))
tplot(p,t);

#Analysis
κymax = 0.0004
steps = 40
σ = 4
N = -σ*A*[1,2,3]
ε = zeros(length(N),steps+1)
κy = zeros(length(N),steps+1)
Nr = zeros(length(N),steps+1)
Myr = zeros(length(N),steps+1)
kax = zeros(length(N),steps+1)
kfl = zeros(length(N),steps+1)
tol = 1e-3
for j in eachindex(N)
    e = [-0.01,0.,0.]
    Δκy = 0.
    for i in 1:steps+1
        s, ks = meso_2d_nw(e,params,Mesh)
        e[1] += (N[j] - s[1] - ks[1,3]*Δκy)/ks[1,1]
        e[3] += Δκy
        error = 1
        iter = 1
        while error > tol
            s, ks = meso_2d_nw(e,params,Mesh)
            error = abs(N[j] - s[1])/abs(N[j])
            e[1] += (N[j] - s[1])/ks[1,1]
            if iter > 20
                println("No convergence")
                break
            end
            iter += 1
        end
        ε[j,i] += e[1]
        κy[j,i] += e[3]
        Nr[j,i] += s[1]
        Myr[j,i] += s[3]
        kax[j,i] += ks[1,1]
        kfl[j,i] += ks[3,3]
        if i == 1; Δκy += κymax/steps; end
    end
end

#Plots
plt.figure(figsize=(8,3))
plt.subplot(1,2,1)
plt.plot(κy[1,:],Myr[1,:]/1e6,"--k",label="N = $(round(Int,N[1]/1e3)) kN")
plt.plot(κy[2,:],Myr[2,:]/1e6,"-k",label="N = $(round(Int,N[2]/1e3)) kN")
plt.plot(κy[3,:],Myr[3,:]/1e6,"-.k",label="N = $(round(Int,N[3]/1e3)) kN")
plt.xlabel(L"κ_y \, [1/mm]")
plt.ylabel(L"M_y \, [kN-m]")
plt.legend()
plt.subplot(1,2,2)
plt.plot(κy[1,:],kfl[1,:]/1e9,"--k",label="N = $(round(Int,N[1]/1e3)) kN")
plt.plot(κy[2,:],kfl[2,:]/1e9,"-k",label="N = $(round(Int,N[2]/1e3)) kN")
plt.plot(κy[3,:],kfl[3,:]/1e9,"-.k",label="N = $(round(Int,N[3]/1e3)) kN")
plt.xlabel(L"κ_y \, [1/mm]")
plt.ylabel(L"\partial M_y/\partial\kappa_y \, [kN-m^2]")
plt.legend();

lstep = steps + 1
l1 = 1
l2 = 2
umax = 32.
umin = 0.01
for step in 1:20:lstep
    ut1, λt1 = meso_constrainedsolver(p,t,en,(x,y)->ffun(x,y;params,ε=ε[l1,step],κx=0,κy=κy[l1,step]),k)
    ut2, λt2 = meso_constrainedsolver(p,t,en,(x,y)->ffun(x,y;params,ε=ε[l2,step],κx=0,κy=κy[l2,step]),k)
#     umax = max(maximum(ut1),maximum(ut2))
#     umin = min(minimum(ut1),minimum(ut2))
    
    fig = plt.figure(figsize=(10,3))
    subfigs = fig.subfigures(1,2,wspace=0.01,width_ratios=[1,2])
    subfigs[1].subplots(1,1)
    plt.plot(κy[l1,1:lstep],Myr[l1,1:lstep]/1e6,"--k",label="P/A = $(round(-N[l1]/A)) MPa")
    plt.plot(κy[l1,step],Myr[l1,step]/1e6,"or")
    plt.plot(κy[l2,1:lstep],Myr[l2,1:lstep]/1e6,"-k",label="P/A = $(round(-N[l2]/A)) MPa")
    plt.plot(κy[l2,step],Myr[l2,step]/1e6,"or")
    plt.gca().ticklabel_format(axis="x",style="sci",scilimits=(0,0))
    plt.xlabel(L"κ_y \, [1/mm]")
    plt.ylabel(L"M_y \, [kN-m]")
    plt.legend()
    axs = subfigs[2].subplots(1,2)
    plt.subplot(1,2,1)
    tplot(p,t,ut1;umin,umax,cbar=false)
    plt.xlabel(L"x \, [mm]")
    plt.ylabel(L"y \, [mm]")
    plt.title("P/A = $(round(Int,-N[l1]/A)) MPa")
    plt.subplot(1,2,2)
    tplot(p,t,ut2;umin,umax,cbar=false)
    plt.gca().get_yaxis().set_visible(false)
    plt.xlabel(L"x \, [mm]")
    plt.title("P/A = $(round(Int,-N[l2]/A)) MPa")
    subfigs[1].subplots_adjust(bottom=0.16,top=0.9,left=0.21,right=0.9,wspace=0.15)
    subfigs[2].subplots_adjust(bottom=0.16,top=0.9,left=0.07,right=1.,wspace=0.15)
    cbar = plt.colorbar(ax=axs[:],ticks=LinRange(umin,umax,5),format="%.0f",
        aspect=25,location="right",pad=0.02,label=L"p \, [MPa]")
#     savefig("mccomp-$(step).pdf")
end;

lstep = steps + 1
l1 = 2
umax = 32.
umin = 0.01
for step in 1:20:lstep
    ut1, λt1 = meso_constrainedsolver(p,t,en,(x,y)->ffun(x,y;params,ε=ε[l1,step],κx=0,κy=κy[l1,step]),k)
#     umax = max(maximum(ut1),maximum(ut2))
#     umin = min(minimum(ut1),minimum(ut2))
    
    fig = plt.figure(figsize=(10,3))
    plt.subplot(1,3,1)
    plt.plot(κy[l1,1:lstep],Myr[l1,1:lstep]/1e6,"-k",label="P/A = $(round(-N[l1]/A)) MPa")
    plt.plot(κy[l1,step],Myr[l1,step]/1e6,"or")
    plt.gca().ticklabel_format(axis="x",style="sci",scilimits=(0,0))
    plt.xlabel(L"κ_y \, [1/mm]")
    plt.ylabel(L"M_y \, [kN-m]")
    plt.legend()
    plt.subplot(1,3,2)
    plt.plot(κy[l1,1:lstep],kfl[l1,1:lstep]/1e9,"-k",label="P/A = $(round(-N[l1]/A)) MPa")
    plt.plot(κy[l1,step],kfl[l1,step]/1e9,"or")
    plt.gca().ticklabel_format(axis="x",style="sci",scilimits=(0,0))
    plt.xlabel(L"κ_y \, [1/mm]")
    plt.ylabel(L"\partial M_y/\partial\kappa_y \, [kN-m^2]")
    plt.legend()
    plt.subplot(1,3,3)
    tplot(p,t,ut1;umin,umax,cbar=false)
    plt.xlabel(L"x \, [mm]")
    plt.ylabel(L"y \, [mm]")
    plt.subplots_adjust(bottom=0.16,top=0.9,left=0.07,right=0.97,wspace=0.3)
    cbar = plt.colorbar(ticks=LinRange(umin,umax,5),format="%.0f",
        aspect=25,location="right",pad=0.05,label=L"p \, [MPa]")
#     savefig("mc-$(step).pdf")
end;


end