export randominside, particle_from_coordinates, coordinates_from_particle

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
Returns the tuple `(p, s, phi)` being `p` of type `Particle` located randomly in the boundary of the unit cell with Birkhoff coordinates `s`, `phi`. If sides is true it also returns the side of the cell where the particle is located
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

"""
Given the birkhof coordinates of a particle `coord` located in the border of the unit cell returns a `Particle` with the corresponding Cartesian coordinates and velocity
"""
function particle_from_coordinates(coord::Vector{T}, n::Int, bt::Vector{<:Obstacle{T}}) where {T <: AbstractFloat}
    R = norm(bt[1].ep - bt[1].sp)
    s_max = R*n   ##It is more efficient to compute this just once
    s_0 = coord[1]
    side = ceil(Int, s_0*n)
    s_0 *= s_max
    real_d = s_0 - (side-1.)*R
    dist = (bt[side].ep - bt[side].sp)
    dist /= norm(dist)
    xp, yp = bt[side].sp + dist*real_d
    
    position = SVector{2, T}([xp, yp])

    sinphi = coord[2]
    phi0 = asin(sinphi)
    vel0 = bt[side].normal/norm(bt[side].normal) ##Take this as the vector t

    vel = SVector{2, T}(rotation(phi0, vel0))
    current_cell = zeros(vel)

    Particle(position, vel, current_cell)

end


"""
Given a particle and the side where it is located, it returns the birkhoff coordinates.
"""
function coordinates_from_particle(p::Particle{T}, n::Int64, bt::Vector{<:Obstacle{T}}, side::Int64) where {T<:AbstractFloat}
    R = norm(bt[1].ep - bt[1].sp) #To do: avoid computing this always
    coord = p.pos + p.current_cell
    s = norm(coord - bt[side].sp) + (side-1)*R
    s /= R*n
    sinphi = DynamicalBilliards.cross2D(bt[side].normal, p.vel)/norm(bt[side].normal)

    return [s, sinphi]

end
    

#    ##Solo lo uso para las condiciones uniformes
# function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int; lyapunov = false, algorithm = randominside)
#     if lyapunov
#             birk_coord = zeros(N,4)
#     else
#         birk_coord = zeros(N,3)
#     end

#     ###initialize
#     p,s,phi = randominside(polygon_sides, bt)

#     birk_coord[1, 1:2] = [s,phi]
#     if lyapunov
#         exps, pos = lyapunovspectrum(p, bt, t; displacement = true)
#         birk_coord[1, 3] = norm(pos[end] -pos[1])
#         birk_coord[1,4] = exps[1]
#     else
#         pos = displacement!(p, bt, t)
#         birk_coord[1, 3] = norm(pos[end] -pos[1])
#     end

    
#     ####
    
#     for i in 2:N
#         p,s,phi = algorithm(polygon_sides, bt)
#         birk_coord[i, 1:2] = [s,phi]
#         if lyapunov
#             exps, pos = lyapunovspectrum(p, bt, t; displacement = true)
#             birk_coord[i, 3] = norm(pos[end] -pos[1])
#             birk_coord[i,4] = exps[1]
#         else
#             pos = displacement!(p, bt, t)
#             birk_coord[i, 3] = norm(pos[end] -pos[1])
#         end
#     end
#     birk_coord
# end


