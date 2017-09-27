import randominside

export randominside, birkhoff_visual, particle_from_coordinates, coordinates_from_particle

"""
Rotates a given vector `vec` by an angle `theta`
"""
function rotation(theta::T, vec::Union{Vector{T}, SVector{2,T}}) where {T <: AbstractFloat}
    q = similar(vec)
    q[1] = cos(theta)*vec[1] - sin(theta)*vec[2]
    q[2] = sin(theta)*vec[1] + cos(theta)*vec[2]
    q
end
"""
Returns the tuple `(p, s, phi)` being `p` of type `Particle` located randomly in the boundary of the unit cell with Birkhoff coordinates `s`, `phi`.
"""
function randominside(bt::Vector{<:Obstacle}, n::Int64; sides = false)
    R = norm(bt[1].ep - bt[1].sp)
    s_max = R*n
    s_0 = rand()*s_max
    side = ceil(Int, s_0*n/s_max)
    real_d = s_0 - (side-1.)*R
    dist = (bt[side].ep - bt[side].sp)
    dist /= norm(dist)
    coord = bt[side].sp + dist*real_d

    sinphi = rand()*2-1.

    phi0 = asin(sinphi)
    vel0 = bt[side].normal/norm(bt[side].normal) ##Take this as the vector to be rotated
    vel = SVector{2, Float64}(rotation(phi0, vel0))

    current_cell = zeros(vel)

    if sides
        return  Particle(coord, vel, current_cell), s_0/s_max, sinphi, side
    end

    Particle(coord, vel, current_cell), s_0/s_max, sinphi

end


function particle_from_coordinates(c::Vector, n::Int, bt::Vector{Obstacle})
    R = norm(bt[1].ep - bt[1].sp)
    s_max = R*n   ##It is more efficient to compute this just once
    s_0 = c[1]
    side = ceil(Int, s_0*n/s_max)
    real_d = s_0 - (side-1.)*R
    dist = (bt[side].ep - bt[side].sp)
    dist /= norm(dist)
    xp, yp = bt[side].sp + dist*real_d

    sinphi = c[2]
    phi0 = asin(sinphi)
    vel0 = bt[side].normal/norm(bt[side].normal) ##Take this as the vector t

    vel = rotation(phi0, vel0)
    current_cell = zeros(vel)

    Particle([xp, yp], vel, current_cell)

end



function coordinates_from_particle(p::Particle, n::Int, bt::Vector{Obstacle}, side::Int)
    R = norm(bt[1].ep - bt[1].sp) #To do: avoid computing this always
    coord = p.pos + p.current_cell
    s = norm(coord - bt[side].sp) + (side-1)*R
    

    sinphi = cross2D(bt[side].normal, p.vel)/norm(bt[side].normal)

    return [s, sinphi]

end

    

   ##Solo lo uso para las condiciones uniformes
function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int; lyapunov = false, algorithm = randominside)
    if lyapunov
            birk_coord = zeros(N,4)
    else
        birk_coord = zeros(N,3)
    end

    ###initialize
    p,s,phi = randominside(polygon_sides, bt)

    birk_coord[1, 1:2] = [s,phi]
    if lyapunov
        exps, pos = lyapunovspectrum(p, bt, t; displacement = true)
        birk_coord[1, 3] = norm(pos[end] -pos[1])
        birk_coord[1,4] = exps[1]
    else
        pos = displacement!(p, bt, t)
        birk_coord[1, 3] = norm(pos[end] -pos[1])
    end

    
    ####
    
    for i in 2:N
        p,s,phi = algorithm(polygon_sides, bt)
        birk_coord[i, 1:2] = [s,phi]
        if lyapunov
            exps, pos = lyapunovspectrum(p, bt, t; displacement = true)
            birk_coord[i, 3] = norm(pos[end] -pos[1])
            birk_coord[i,4] = exps[1]
        else
            pos = displacement!(p, bt, t)
            birk_coord[i, 3] = norm(pos[end] -pos[1])
        end
    end
    birk_coord
end


