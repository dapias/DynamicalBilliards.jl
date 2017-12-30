export lyapunovmaximum, distance, xdistance

import Base.copy 
function copy(p::Particle{T}) where {T<:AbstractFloat}
    pos = copy(p.pos)
    vel = copy(p.vel)
    cc = copy(p.current_cell)

    Particle(pos,vel, cc)
end
    

function lyapunovmaximum(particle::Particle{T}, bt::Vector{Obstacle{T}}, t::T) where {T<:AbstractFloat}
    p = copy(particle)
    lyapmax = lyapunovspectrum!(p,bt, t)[1]
end

function distance(particle::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T) where {T<:AbstractFloat}
    p = copy(particle)
    rpos = SVector{2,T}[]
    push!(rpos, p.pos)
    count = zero(T)

    while count < t
        tmin::T, i::Int = next_collision(p, bt)
        # set counter
        if count +  increment_counter(t, tmin) > t
            break
        end
        #####
        tmin = relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        count += increment_counter(t, tmin)
    end#time loop

    tmin = t - count 
    propagate!(p, tmin)
    push!(rpos, p.pos + p.current_cell)
    d = norm(rpos[2]-rpos[1])
end

function xdistance(particle::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T) where {T<:AbstractFloat}
    p = copy(particle)
    rpos = SVector{2,T}[]
    push!(rpos, p.pos)
    count = zero(T)

    while count < t
        tmin::T, i::Int = next_collision(p, bt)
        # set counter
        if count +  increment_counter(t, tmin) > t
            break
        end
        #####
        tmin = relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        count += DynamicalBilliards.increment_counter(t, tmin)
    end#time loop

    tmin = t - count 
    propagate!(p, tmin)
    push!(rpos, p.pos + p.current_cell)
    #    d = norm(rpos[2]-rpos[1])
    d = rpos[2][1] - rpos[1][1] #distance in x
end
