export randominside, birkhoff_visual

function rotation(theta::Float64, vec::SVector{2,Float64})
    q = similar(vec)
    q[1] = cos(theta)*vec[1] - sin(theta)*vec[2]
    q[2] = sin(theta)*vec[1] + cos(theta)*vec[2]
    q
end
"""
Returns the tuple `(p, s, phi)` being `p` of type `Particle` located randomly in the boundary of the unit cell with Birkhoff coordinates `s`, `phi`.
"""
function randominside(n::Int, bt::Vector{Obstacle})
    R = norm(bt[1].ep - bt[1].sp)
    s = R*n
    initial_point = rand()*s
    side = ceil(Int, initial_point*n/s)
    real_d = initial_point - (side-1.)*R
    dist = (bt[side].ep - bt[side].sp)
    dist /= norm(dist)
    xp, yp = bt[side].sp + dist*real_d

    sinphi = rand()*2-1.
    phi0 = asin(sinphi)
    vel0 = bt[side].normal/norm(bt[side].normal) ##Take this as the vector to be rotated
    vel = rotation(phi0, vel0)

    current_cell = zeros(vel)

    Particle([xp, yp], vel, current_cell), initial_point, sinphi

end

function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int; lyapunov = false)
    if lyapunov
            birk_coord = zeros(N,4)
    else
        birk_coord = zeros(N,3)
    end
    
    for i in 1:N
        p,s,phi = randominside(polygon_sides, bt)
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
