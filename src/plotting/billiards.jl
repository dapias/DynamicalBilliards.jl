using PyPlot, StaticArrays
export plot_billiard, billiard_julia, plot_triangular


"""
```julia
plot_billiard(bt::Vector{Obstacle})
```
Plot all obstacles in `bt` using the default arguments, set
`xlim` and `ylim` to be 10% larger than `cellsize` and
set the axis aspect ratio to equal.

```julia
plot_billiard(bt, xmin, ymin, xmax, ymax)
```
Plot the given (periodic) billiard `bt` on the current PyPlot figure, repeatedly
plotting from `(xmin, ymin)` to `(xmax, ymax)`. Only works for rectangular billiards.

```julia
plot_billiard(bt, xt::Vector, yt::Vector; plot_orbit = true)
```
Plot the given (periodic) billiard `bt` along with a particle trajectory defined
by `xt` and `yt`, on the current PyPlot figure. Only works for rectangular billiards.

Sets limits automatically. Set the keyword argument `plot_orbit = false` to not
plot the orbit defined by `(xt, yt)`.
"""
function plot_billiard(bt::Vector{<:Obstacle{T}}) where {T}
  for obst in bt
    plot_obstacle(obst)
  end
  xmin, ymin, xmax, ymax = cellsize(bt)
  dx = xmax - xmin; dy = ymax - ymin
  PyPlot.xlim(xmin - 0.1dx, xmax + 0.1dx)
  PyPlot.ylim(ymin - 0.1dy, ymax + 0.1dy)
  PyPlot.gca()[:set_aspect]("equal")
end


function plot_billiard(bt, xmin, ymin, xmax, ymax)
  # Cell limits:
  cellxmin, cellymin, cellxmax, cellymax = cellsize(bt)
  dcx = cellxmax - cellxmin
  dcy = cellymax - cellymin
  # Obstacles to plot:
  toplot = Obstacle{eltype(bt)}[]
  for obst in bt
    typeof(obst) <: PeriodicWall && continue
    push!(toplot, obst)
  end
  # Find displacement vectors (they will multiply dcx, dcy)
  dx = (floor((xmin - cellxmin)/dcx):1:ceil((xmax - cellxmax)/dcx))*dcx
  dy = (floor((ymin - cellymin)/dcy):1:ceil((ymax - cellymax)/dcy))*dcy
  # Plot displaced Obstacles
  for x in dx
    for y in dy
      disp = SVector(x,y)
      for obst in toplot
        plot_obstacle(translation(obst, disp))
      end
    end
  end
  # Set limits etc.
  PyPlot.xlim(xmin, xmax)
  PyPlot.ylim(ymin, ymax)
  PyPlot.gca()[:set_aspect]("equal")
end

function plot_billiard(bt, xt::AbstractVector, yt::AbstractVector; plot_orbit = true)
  xmin = floor(minimum(round.(xt,3))); xmax = ceil(maximum(round.(xt,3)))
  ymin = floor(minimum(round.(yt,3))); ymax = ceil(maximum(round.(yt,3)))
  if plot_orbit
    plot(xt, yt, color = "blue")
  end

  plot_billiard(bt, xmin, ymin, xmax, ymax)
end


"""
    translation(obst::Obstacle, vector)
Create a copy of the given obstacle with its position
translated by `vector`.
"""
function translation(d::Circular, vec)
  newd = typeof(d)(d.c .+ vec, d.r)
end

function translation(w::Wall, vec)
  neww = typeof(w)(w.sp + vec, w.ep + vec, w.normal)
end


"""
```julia
billiard_julia(; plotit = true)
```
Return the awesome "Julia-logo" billiard shown in the introduction
of DynamicalBilliards.jl.

By default it also plots the billiard in a new `PyPlot.figure()` using the correct colors.
"""
function billiard_julia(; plotit = true)

  bt = billiard_rectangle()

  r = 0.165
  ewidth = 6.0
  redcent = [0.28, 0.32]
  red = Disk(redcent, r, "Red dot")
  purple = Disk([1 - redcent[1], redcent[2]], r, "Purple dot")
  green = Disk([0.5, 1 - redcent[2]], r, "Green dot")
  push!(bt, red, purple, green)

  if plotit == true
    PyPlot.figure()
    for w in bt
      plot_obstacle(w; color = (0,0,0, 1), linewidth = 3.0)
    end
    plot_obstacle(red; edgecolor = (203/255, 60/255, 51/255),
    facecolor = (213/255, 99/255, 92/255), linewidth = ewidth)
    plot_obstacle(purple; edgecolor = (149/255, 88/255, 178/255),
    facecolor = (170/255, 121/255, 193/255), linewidth = ewidth)
    plot_obstacle(green, edgecolor = (56/255, 152/255, 38/255),
    facecolor = (96/255, 173/255, 81/255), linewidth = ewidth)

    # particle edge color
    # darkblue = (64/255, 99/255, 216/255)
    # lightblue = (102/255, 130/255, 223/255)

    PyPlot.axis("off")
    PyPlot.tight_layout()
    PyPlot.gca()[:set_aspect]("equal")
    PyPlot.xlim(-0.1,1.1)
    PyPlot.ylim(-0.1,1.1)
  end

  return bt
end

"""
Draw a triangular periodic billiard with a hexagon as unit cell. The variable `space` passed is the space betweem the scatterrers
"""

function plot_triangular(bt, xmin, ymin, xmax, ymax, space) 
    plot_billiard(bt)
    disp1 = [0.0, space]
    disp2 = [space*cos(pi/6), space*sin(pi./6)]
    
    PyPlot.xlim(xmin , xmax)
    PyPlot.ylim(ymin , ymax)
    PyPlot.gca()[:set_aspect]("equal")
    
    d = bt[find(x -> isa(x,Circular), bt)][1]
    
    j =0
    while xmin < -j*space*cos(pi/6) + d.c[1]
        d_temp = translation(d, -j*disp2)
        plot_obstacle(d_temp)
        i = 1
        while ymin < -i*space + d_temp.c[2]
            plot_obstacle(translation(d_temp, -i*disp1))
            i +=1 
        end
        
        k = 1
        while ymax > k*space + d_temp.c[2]
            plot_obstacle(translation(d_temp, k*disp1))
            k +=1 
        end
        j +=1 
    end
    
    
    j = 1
     while xmax > j*space*cos(pi/6)+ d.c[1]
        d_temp = translation(d, j*disp2)
        plot_obstacle(d_temp)
        i = 1
        while ymin < -i*space + d_temp.c[2]
            plot_obstacle(translation(d_temp, -i*disp1))
            i +=1 
        end
        
        k = 1
        while ymax > k*space + d_temp.c[2]
            plot_obstacle(translation(d_temp, k*disp1))
            k +=1 
        end
            j +=1 
      end
end

function plot_triangular(bt, xt::Vector{Float64}, yt::Vector{Float64}; plot_orbit = true)
  
  space = norm(bt[1].sp - bt[1].ep)*cos(pi/6)*2
    
  xmin = floor(minimum(round.(xt,8)))  - space
  xmax = ceil(maximum(round.(xt,8))) + space
  ymin = floor(minimum(round.(yt,8))) - space
  ymax = ceil(maximum(round.(yt,8))) +  space
    
  if plot_orbit
    plot(xt, yt, color = "blue")
  end
    plot_triangular(bt, xmin, ymin, xmax, ymax, space)
end
