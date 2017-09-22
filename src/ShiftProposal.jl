export birkhoff_visual, motion_MCMC, acceptance, shift_evolution

mutable struct ShiftProposal <: Proposal
    parameter::Union{Function, Float64}    #t_shift
    f::Function
end


function shift_evolution(particle::Particle, bt::Vector{Obstacle}, time::Float64, n::Int64, index::Int64)

    p = copy(particle)
    if time < 0.0
        p.vel = -p.vel
        p.pos += bt[index].normal 
    end

    t = abs(time)
    count = 0.0
    colobst_idx = 1
    t_to_write = 0.0
    
    while count < t
        tcol::Float64 = 0.0
        tmin::Float64 = Inf
        for i in eachindex(bt)
            tcol = collisiontime(p, bt[i])
            # Set minimum time:
            if tcol < tmin
                tmin = tcol
                colobst_idx = i
            end
        end
        
        if count + increment_counter(t,tmin) > t
            break
        else
            count += increment_counter(t, tmin)
        end
        
        propagate!(p, tmin)
        resolvecollision!(p, bt[colobst_idx])
        t_to_write += tmin
        if typeof(bt[colobst_idx]) <: PeriodicWall
            continue
        else
            t_to_write = 0.0
        end
    end#time loop
    
    tmin = t - count 
    propagate!(p, tmin)
    ##Here the evolution ends and then we propagate until the next collision
    new_index = evolution_to_periodicwall!(p, bt, n)
    new_particle = Particle(p.pos, p.vel, [0.,0.])


    if time < 0.0
        new_particle.vel = -new_particle.vel
        new_particle.pos += bt[new_index].normal
        new_index = new_index >= 4 ? mod(new_index + div(n,2), n) : new_index + div(n,2)

    end

    return new_particle, new_index
    
end

"""
Evolves the particle until it reaches the first periodic wall
"""
function evolution_to_periodicwall!(p::Particle, bt::Vector{Obstacle}, n::Int64)

    colobst_idx = 1
    
    while true
        tcol::Float64 = 0.0
        tmin::Float64 = Inf

        for i in eachindex(bt)
            tcol = collisiontime(p, bt[i])
            if tcol < tmin
                tmin = tcol
                colobst_idx = i
            end
        end#obstacle loop
        
        propagate!(p, tmin)
        resolvecollision!(p, bt[colobst_idx])
        
        if typeof(bt[colobst_idx]) <: PeriodicWall
            break
        end

    end

    index = colobst_idx >= 4 ? mod(colobst_idx + div(n,2), n) : colobst_idx + div(n,2)


    return index
end
    
function proposal(init_particle::Particle, sides::Int64, bt::Vector{Obstacle}, tshift::Float64, index::Int64)

    particle, ind =  shift_evolution(init_particle, bt, tshift, sides, index)
    
end

#beta =0
function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64, tshift::Float64; prop = "local")

    pos = zeros(N,2)
    ###initialize
    p,s, sinphi, index = randominside(polygon_sides, bt, sides = true)
    pos[1,:] = [s, sinphi]
    ####
    ##Needed types
    proposal_dist1 = ShiftProposal(tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist1.parameter, i) )
    proposal_dist2 = ShiftProposal(-tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist2.parameter, i) )
    ###############
    for i in 2:N
        proposal_dist = proposal_dist1
        if rand() < 0.5
            proposal_dist = proposal_dist2
        end
        
        k, ind =  proposal_dist.f(p, index)
        s, sinphi = coordinates_from_particle(k, polygon_sides, bt, ind)
        pos[i, :] = [s, sinphi]
        p, index = k, ind
    end
    pos
end


function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64, tshift::Float64, prop ::String)

    birk_coord = zeros(N,3)
    ###initialize
    p,s, sinphi, index = randominside(polygon_sides, bt, sides = true)
    birk_coord[1, 1:2] = [s,sinphi]
    dist = distance(p, bt, t)
    birk_coord[1, 3] = dist    
    ####
    ##Needed types
    proposal_dist1 = ShiftProposal(tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist1.parameter, i) )

    proposal_dist2 = ShiftProposal(-tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist2.parameter, i) )
    ###############
    for i in 2:N
        proposal_dist = proposal_dist1
        random = rand()
        if rand() < 0.5
            proposal_dist = proposal_dist2
        end
        k, ind =  proposal_dist.f(p, index)
        s, sinphi = coordinates_from_particle(k, polygon_sides, bt, ind)
        dist = distance(p, bt, t)
        birk_coord[i, :] = [s,sinphi, dist]
        p, index = k, ind
    end
    
    birk_coord
end

#beta no 0
function acceptance(p::Particle, s::Float64, sinphi::Float64, polygon_sides::Int, bt::Vector{Obstacle}, prop_dist::Proposal, cano::CanonicalDistribution,  time::Float64, index::Int64)
    
    E1 = cano.f(p)
    p1, ind = prop_dist.f(p, index)
    E2 = cano.f(p1)
    a = 0.0
    if E2 - E1 >= 0.0
        p = p1
        E1 = E2
        index = ind
    else
        a = exp(E2 - E1)
        r = rand()
        if r < a   #If r is less than a the motion is accepted
            p = p1
            E1 = E2
            index = ind
        end
    end

    dist = -E1/cano.beta
    s, sinphi = coordinates_from_particle(p, polygon_sides, bt, index)
    
    return p, s, sinphi, dist, index

end


function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64, beta::Float64, tshift::Float64, time::Float64, prop::String)
    ##Needed types
    proposal_dist1 = ShiftProposal(tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist1.parameter, i) )
    proposal_dist2 = ShiftProposal(-tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist2.parameter, i) )
    obs = x::Particle -> distance(x, bt, time)
    canonical = CanonicalDistribution(beta, obs)
    #######
    pos = zeros(N,2)
    p,s, sinphi, index = randominside(polygon_sides, bt, sides = true)
    pos[1,:] = [s, sinphi]

    for i in 2:N
        proposal_dist = proposal_dist1
        random = rand()
        if rand() < 0.5
            proposal_dist = proposal_dist2
        end
        p, s, sinphi, dist, index = acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, canonical,  time, index)
        pos[i, :] = [s, sinphi]
    end
    pos
end

 function birkhoff_visual(time::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64, beta::Float64, tshift::Float64, prop::String)

     birk_coord = zeros(N,3)
#     ###initialize
     p,s, sinphi, index = randominside(polygon_sides, bt, sides = true)
     birk_coord[1, 1:2] = [s,sinphi]
     dist = distance(p, bt, time)
     birk_coord[1, 3] = dist    
#     ####

    ##Needed types
    proposal_dist1 = ShiftProposal(tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist1.parameter, i) )
     proposal_dist2 = ShiftProposal(-tshift, (x::Particle, i::Int64) -> proposal(x, polygon_sides, bt, proposal_dist2.parameter, i) )
    obs = x::Particle -> distance(x, bt, time)
     canonical = CanonicalDistribution(beta, obs)
    #######
    
     for i in 2:N
         proposal_dist = proposal_dist1
         random = rand()
         if rand() < 0.5
            proposal_dist = proposal_dist2
         end
         p, s, sinphi, dist, index = acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, canonical,  time, index)
         birk_coord[i, :] = [s,sinphi, dist]
     end
    
     birk_coord
end
