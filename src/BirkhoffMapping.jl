export randominside, particle_from_coordinates, coordinates_from_particle, InitialCondition


"""
Type that contains the particle, its Birkhoff coordinates and an index representing the side of the polygon where it is located
"""
struct InitialCondition{T<:AbstractFloat}
    particle::Particle{T}
    s::T
    sinphi::T
    index::Int
end
"""
Rotates a given vector `vec` by an angle `theta`
"""
function rotation(theta::T, vec::Union{Vector{T}, SVector{2,T}}) where {T <: AbstractFloat}
    #q = similar(vec)
    q = Vector{T}(2)
    q[1] = cos(theta)*vec[1] - sin(theta)*vec[2]
    q[2] = sin(theta)*vec[1] + cos(theta)*vec[2]
    q
end
"""
Returns the tuple `(p, s, phi)` being `p` of type `Particle` located randomly in the boundary of the unit cell with Birkhoff coordinates `s`, `phi`. If sides is true it also returns the side of the cell where the particle is located
"""
function randominside(bt::Vector{<:Obstacle{T}}, n::Int) where {T<:AbstractFloat}
    R = norm(bt[1].ep - bt[1].sp)
    s_max = R*n
    s_0 = rand()*s_max
    side = ceil(Int, s_0*n/s_max)
    real_d = s_0 - (side-1.)*R
    dist = (bt[side].ep - bt[side].sp)
    dist /= norm(dist)
    coord = bt[side].sp + dist*real_d

    sinphi = T(rand()*2-1.)

    phi0 = asin(sinphi)
    vel0 = T.(bt[side].normal/norm(bt[side].normal)) ##Take this as the vector to be rotated
    vel = SVector{2, T}(rotation(phi0, vel0))
    current_cell = zeros(vel)
    
    return  InitialCondition(Particle(coord, vel, current_cell), s_0/s_max, sinphi, side)    
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
function coordinates_from_particle(p::Particle{T}, n::Int, bt::Vector{<:Obstacle{T}}, side::Int) where {T<:AbstractFloat}
    R = norm(bt[1].ep - bt[1].sp) #To do: avoid computing this always
    coord = p.pos + p.current_cell
    s = norm(coord - bt[side].sp) + (side-1)*R
    s /= R*n
    sinphi = DynamicalBilliards.cross2D(bt[side].normal, p.vel)/norm(bt[side].normal)

    return [s, sinphi]

end
