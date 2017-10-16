# ParticlesObstacles.jl must be loaded BEFORE this
export lyapunovspectrum!

##Auxiliar Functions ##
"""
```julia
gramschmidt(u::Matrix{T})
```
Apply the Gram-Schmidt procedure to a square matrix of 4  vectors 
"""
function gramschmidt(u::Matrix{T}) where {T}
    w = eye(T, 4)
    w[:,1] = u[:,1];
    v1 = w[:,1]/norm(w[:,1])
    w[:,2] = u[:,2] - dot(u[:,2],v1)*v1;
    v2 = w[:,2]/norm(w[:,2]);
    w[:,3] = (u[:,3] - dot(u[:,3],v2)*v2 - dot(u[:,3],v1)*v1)
    v3 = w[:,3]/norm(w[:,3])
    w[:,4] = u[:,4] - dot(u[:,4],v3)*v3 - dot(u[:,4],v2)*v2 - dot(u[:,4],v1)*v1 
    return w
end
"""
    specular!(p::AbstractParticle, o::Obstacle, offset::Matrix)
Perform specular reflection based on the normal vector of the Obstacle.
The function updates the position and velocity of the particle
together with the components of 4 offset vectors stored in the matrix
`offset` as columns.
"""
function specular!(p::AbstractParticle, o::Disk, offset::Matrix)
    n = normalvec(o, p.pos)
    ti = [-p.vel[2],p.vel[1]]
    cosa = dot(n, -p.vel)
    p.vel = p.vel - 2*dot(n, p.vel)*n
    tf = [-p.vel[2], p.vel[1]]
    for k in 1:4
        x = offset[:,k]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 27)
        # with norm(p) = 1
        x[3:4]  = x[3:4] - 2.*dot(x[3:4],n)*n-2./o.r*dot(x[1:2],ti)/cosa*tf
        x[1:2]  = x[1:2] - 2.*dot(x[1:2],n)*n
        ###
        offset[:,k] = x
    end
end


function specular!(p::AbstractParticle, o::FiniteWall, offset::Matrix)
    n = normalvec(o, p.pos)
    specular!(p, o)
    for k in 1:4
        x = offset[:,k]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 20)
        x[1:2]  = x[1:2] -  2.*dot(x[1:2],n)*n
        x[3:4]  = x[3:4] - 2.*dot(x[3:4],n)*n
        ###
        offset[:,k] = x
    end
end
"""
    resolvecollision!(p::AbstractParticle, o::Union{Disk, FiniteWall}, offset::Matrix)
Resolve the collision between particle `p` and obstacle `o` of type *Circular*,
updating the components of the offset vectors stored in the matrix `offset` as columns.
"""
function resolvecollision!(p::AbstractParticle, o::Union{Disk, FiniteWall}, offset::Matrix)::Void
    specular!(p, o, offset)
    return
end

resolvecollision!(p::AbstractParticle, o::PeriodicWall, offset::Matrix)::Void =
resolvecollision!(p, o)

"""
```julia
propagate!(p::AbstractParticle, t, offset::Matrix)
```
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields together with the components of the offset vectors stored in the `offset` matrix.
"""
function propagate!(p::Particle, t::Real, offset::Matrix)
    vx0=p.vel[1]
    vy0=p.vel[2]
    p.pos += [vx0*t, vy0*t]
    for k in 1:4
        x = offset[:,k]
        temp = [x[1] + t*x[3], x[2] + t*x[4], x[3], x[4]]
        offset[:,k] = temp
    end
end

"""
```julia
    relocate(p::AbstractParticle, o::Obstacle, t, offset::Matrix) -> newt
```
Propagate the particle's position for time `t` (corrected) and update the components of the `offset` matrix
"""
function relocate!(p::Particle{T}, o::Obstacle{T}, tmin::T, offset::Matrix{T})::T where {T}
    tmin = relocate!(p, o, tmin)
    for k in 1:4
        x = offset[:,k]
        temp = [x[1] + tmin*x[3], x[2] + tmin*x[4], x[3], x[4]]
        offset[:,k] = temp
    end
    return tmin
end


"""
```julia
    lyapunovspectrum(p::AbstractParticle, bt::Vector{Obstacle}, t::Float64)
```
Returns the finite time lyapunov exponents for a given initial condition of the particle `p` . The time `t` is asked to be of type Float64 .
"""
function lyapunovspectrum!(p::Particle{T}, bt::Vector{Obstacle{T}}, t::T) where {T<:AbstractFloat}
    offset = eye(T, 4) #The unit vectors in the 4 directions

    if t <= 0
        error("`evolve!()` cannot evolve backwards in time.")
    end

    count = zero(T)
    norms = ones(T, 1,4)#Where the norms of the offset vectors will be stored

    while count < t
        # Declare these because `bt` is of un-stable type!
        tmin::T, i::Int = next_collision(p, bt)
        
        # set counter
        if count +  increment_counter(t, tmin) > t
            break
        end
        ###
        tmin = relocate!(p, bt[i], tmin, offset)
        resolvecollision!(p, bt[i], offset)
        count += increment_counter(t, tmin)

        if typeof(bt[i]) <: Wall
            continue
        else
            offset = gramschmidt(offset)
            instantaneous_norms = [norm(offset[:,j]) for j in 1:4]
            norms = vcat(norms,instantaneous_norms')
            for j in 1:4
                offset[:,j] = offset[:,j]/norm(offset[:,j])
            end
        end
    end#time loop

    tmin = t - count
    propagate!(p, tmin, offset)
    instantaneous_norms = [norm(offset[:,j]) for j in 1:4]
    norms = vcat(norms,instantaneous_norms')
    for j in 1:4
        offset[:,j] = offset[:,j]/norm(offset[:,j])
    end

    exps = zeros(T, 4)

    for k in 1:4
        exps[k] = sum(log.(norms[:,k]))/t
   end

    return exps
end
