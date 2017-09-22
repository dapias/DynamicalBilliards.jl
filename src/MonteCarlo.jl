export NeighborhoodProposal, proposal, birkhoff_visual, motion_MCMC, CanonicalDistribution, distance

abstract type Distribution end  #Distributions will be a type
abstract type Proposal end


mutable struct CanonicalDistribution <: Distribution
    beta::Float64
    observable::Function
    f::Function
end

CanonicalDistribution(beta, observable) = CanonicalDistribution(beta, observable,x::Particle -> -beta*observable(x))  

mutable struct NeighborhoodProposal <: Proposal
    parameter::Union{Function, Float64}    #Standard_Deviation in the case of the Gaussian
    f::Function
end



"""
Propose a new point (s, sinphi) based on a half-gaussian distribution g(xprime,x). It rejects if sintheta_new > 1 or < -1
"""
function proposal(init_cond::Vector, sides::Int64, bt::Vector{Obstacle}, sigma::Float64)

    delta_theta = rand()*(2pi)
    rprime = abs(randn()*sigma)
    s_max = sides*norm(bt[1].ep - bt[1].sp)
    delta_s = mod(rprime,s_max)
    s_new = cos(delta_theta)*rprime + init_cond[1]
    s_new = mod(s_new, s_max)
    sintheta_new = sin(delta_theta)*rprime + init_cond[2]

    if abs(sintheta_new) > 1.0
        return init_cond
    end
    
    return [s_new, sintheta_new]
end


"""
Returns the distance that a random particle travels on a billiard table `bt` during the time t
"""
function distance(particle::Particle, bt::Vector{Obstacle}, t::Float64)
    #To do: Valid for a negative time
    p = copy(particle)
    rpos = SVector{2,Float64}[]   #To do: Preinitialize this
    push!(rpos, p.pos)

    count = zero(t)
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
        end#obstacle loop

        # set counter
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
    push!(rpos, p.pos + p.current_cell)

    minus_d = norm(rpos[2]-rpos[1])
    
end


##Functions designed for beta = 0


function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64)

    prop_dist = NeighborhoodProposal(sigma, x::Vector -> proposal(x, polygon_sides, bt, prop_dist.parameter) )
    
    pos = zeros(N,2)
    p,s, sinphi = randominside(polygon_sides, bt)
    pos[1,:] = [s, sinphi]
    for i in 1:N
        s, sinphi = prop_dist.f([s,sinphi])
        pos[i, :] = [s, sinphi]
    end
    pos
end

function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64 )
    
    prop_dist = NeighborhoodProposal(sigma, x::Vector -> proposal(x, polygon_sides, bt, prop_dist.parameter) )
    birk_coord = zeros(N,3)
    
    ###initialize
    p,s,phi = randominside(polygon_sides, bt)
    birk_coord[1, 1:2] = [s,phi]
    birk_coord[1, 3] =  distance(p, bt, t)
       
    ####
    for i in 2:N
        s,phi = prop_dist.f([s,phi])
        birk_coord[i, 1:2] = [s,phi]
        p = particle_from_coordinates([s,phi], polygon_sides, bt)
        birk_coord[i, 3] =  distance(p, bt, t)
    end
    
    birk_coord
end

## Functions for \beta different from zero
function acceptance(p::Particle, s::Float64, sinphi::Float64, polygon_sides::Int, bt::Vector{Obstacle}, prop_dist::Proposal, cano::CanonicalDistribution,  time::Float64)
    
    E1 = cano.f(p)
    s1, sinphi1 = prop_dist.f([s,sinphi])
    p1 = particle_from_coordinates([s1,sinphi1], polygon_sides, bt)
    E2 = cano.f(p1)
    a = 0.0
    if E2 - E1 >= 0.0
        s = s1
        sinphi = sinphi1
        p = p1
        E1 = E2
    else
        a = exp(E2 - E1)
        r = rand()
        if r < a   #If r is less than a the motion is accepted
            s = s1
            sinphi = sinphi1
            p = p1
            E1 = E2
        end
    end

    dist = -E1/cano.beta
    return p, s, sinphi, dist

end

function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64, beta::Float64, time::Float64)

    ##Needed types
    prop_dist = NeighborhoodProposal(sigma, x::Vector -> proposal(x, polygon_sides, bt, prop_dist.parameter) )
    obs = x::Particle -> distance(x, bt, time)
    canonical = CanonicalDistribution(beta, obs)
    ###
    
    pos = zeros(N,2)
    p,s, sinphi = randominside(polygon_sides, bt)
    pos[1,:] = [s, sinphi]
    for i in 1:N
        p,s, sinphi, r = acceptance(p, s, sinphi, polygon_sides, bt, prop_dist, canonical, time)
        pos[i, :] = [s, sinphi]
    end
    pos
end


function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, sigma::Float64, beta::Float64)

    birk_coord = zeros(N,3)
    ###initialize
    p,s, sinphi = randominside(polygon_sides, bt)
    birk_coord[1, 1:2] = [s,sinphi]
    dist = distance(p, bt, t)
    birk_coord[1, 3] = dist    
    ####

    ##Needed types
    proposal_dist = NeighborhoodProposal(sigma, x::Vector -> proposal(x, polygon_sides, bt, proposal_dist.parameter) )
    obs = x::Particle -> distance(x, bt, t)
    canonical = CanonicalDistribution(beta, obs)
    ###############
    
    for i in 2:N
        p, s, sinphi, dist =  acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, canonical, t)
        birk_coord[i, :] = [s,sinphi, dist]
    end
    
    birk_coord
end
