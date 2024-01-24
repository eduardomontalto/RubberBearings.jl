# ---------- PLOTTING FUNCTION ----------
#=
INPUT:
p: array of (x,y) coordinates. Size: nn x 2
t: array of connectivity. Size: nel x 3
u: nodal values of function to be plotted. Size: nn x 1
cmap: colormap ("Spectral_r", "YlOrRd", "YlGnBu_r")
=#
function tplot(p::AbstractArray{Float64}, 
    t::AbstractArray{Int}, 
    u::Union{Nothing,AbstractArray}=nothing;
    cmap::AbstractString="YlOrRd",
    cbar::Bool=true,
    umax::Union{Nothing,Real}=nothing,
    umin::Union{Nothing,Real}=nothing)
plt.axis("equal")
if isnothing(u)
    plt.tripcolor(p[:,1], p[:,2], Array(t .- 1), 0*t[:,1],
              cmap="Set3", edgecolors="k", linewidth=1)
else
    if isnothing(umax); umax = maximum(u); end
    if isnothing(umin); umin = minimum(u); end
    levels = LinRange(umin,umax,50)
    ticks = LinRange(umin,umax,5)
    plt.tricontourf(p[:,1], p[:,2], Array(t .- 1), u[:], levels, cmap=cmap)
    if cbar
        if umax > 100
            plt.colorbar(location="bottom",format="%.0f",ticks=ticks)
        elseif umax > 10
            plt.colorbar(location="bottom",format="%.1f",ticks=ticks)
        elseif umax <= 10
            plt.colorbar(location="bottom",format="%.2f",ticks=ticks)
    end
    end
end
plt.draw()
end;
