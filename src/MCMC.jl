export NeighborhoodProposal, proposal, birkhoff_visual, motion_MCMC, CanonicalDistribution, minus_distance!, lyap_dist

import Distributions.Rayleigh

abstract type Distribution end  #Distributions will be a type
abstract type Proposal end
#abstract type Observable end


mutable struct CanonicalDistribution <: Distribution
    beta::Float64
    observable::Function
    f::Function
end

CanonicalDistribution(beta, observable) = CanonicalDistribution(beta, observable,x::Particle -> exp(-beta*observable(x)))  ##1-1 for extracting the maximum ftle

mutable struct NeighborhoodProposal <: Proposal
    parameter::Union{Function, Float64}    #Standard_Deviation in the case of the Gaussian
    f::Function
end


"""
Propose a new point (s, sinphi) based on a half-gaussian distribution g(xprime,x). It rejects if sintheta_new > 1 or < -1
"""
function proposal(init_cond::Vector, n::Int64, bt::Vector{Obstacle}, sigma::Float64)
    delta_theta = rand()*(2pi)
    rprime = abs(randn()*sigma)
    s_max = n*norm(bt[1].ep - bt[1].sp)
    delta_s = mod(rprime,s_max)
    s_new = cos(delta_theta)*rprime + init_cond[1]
    s_new = mod(s_new, s_max)
    sintheta_new = sin(delta_theta)*rprime + init_cond[2]

    if abs(sintheta_new) > 1.0
        return init_cond
    end
    
    return [s_new, sintheta_new]
end

 

##Function designed for beta = 0
function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, proposal_dist::Proposal; lyapunov = false)
                         
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
        s,phi = proposal_dist.f([s,phi])
        p = particle_from_coordinates([s,phi], polygon_sides, bt)
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


function acceptance(p::Particle, s::Float64, sinphi::Float64, polygon_sides::Int, bt::Vector{Obstacle}, proposal_dist::Proposal, cano::CanonicalDistribution)
    E1 = cano.f(p)
    s1, sinphi1 = proposal_dist.f([s,sinphi])
    p1 = particle_from_coordinates([s1,sinphi1], polygon_sides, bt)
    E2 = cano.f(p1)
    a = min(1.0, E2/E1)
    r = rand()
    if r < a
        s = s1
        sinphi = sinphi1
        p = p1
    end

    return p, s, sinphi

end


##Function designed for beta different from 0
function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, proposal_dist::Proposal, cano::CanonicalDistribution; lyapunov = false)
                         
    if lyapunov
            birk_coord = zeros(N,4)
    else
        birk_coord = zeros(N,3)
    end

    ###initialize
    p,s, sinphi = randominside(polygon_sides, bt)

    birk_coord[1, 1:2] = [s,sinphi]
    if lyapunov
        exps, pos = lyapunovspectrum(p, bt, t; displacement = true)
        birk_coord[1, 3] = norm(pos[end] -pos[1])
        birk_coord[1,4] = exps[1]
    else
        pos = displacement!(p, bt, t)
        birk_coord[1, 3] = norm(pos[end] -pos[1])
    end

    
    ####
    ##Check this part because displacement! changes the particule but lyapunov spectrum not. Then this is a problem. We have to be consistent.
    for i in 2:N
        p, s, sinphi =  acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, cano)
        birk_coord[i, 1:2] = [s,sinphi]
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




function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, proposal_dist::NeighborhoodProposal)
    pos = zeros(N,2)
    p,s, sinphi = randominside(polygon_sides, bt)
    pos[1,:] = [s, sinphi]
    for i in 1:N
        s, sinphi = proposal_dist.f([s,sinphi])
        pos[i, :] = [s, sinphi]
    end
    pos
end
        
function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, proposal_dist::NeighborhoodProposal, cano::CanonicalDistribution)
    pos = zeros(N,2)
    p,s, sinphi = randominside(polygon_sides, bt)
    pos[1,:] = [s, sinphi]
    for i in 1:N
        p,s, sinphi = acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, cano)
        pos[i, :] = [s, sinphi]
    end
    pos
end

##For  the pair(l, d) as observable
function motion_MCMC(N::Int, bt::Vector{Obstacle}, polygon_sides::Int, proposal_dist::NeighborhoodProposal, beta::Float64, time::Float64)
    pos = zeros(N,2)
    p,s, sinphi = randominside(polygon_sides, bt)
    pos[1,:] = [s, sinphi]
    for i in 1:N
        p,s, sinphi, r, lyap = acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, beta, time)
        pos[i, :] = [s, sinphi]
    end
    pos
end



"""
Returns the initial and the final position in an array
"""
function minus_distance!(p::Particle, bt::Vector{Obstacle}, t)
    if t <= 0
        error("`evolve!()` cannot evolve backwards in time.")
    end
    rpos = SVector{2,Float64}[]
    push!(rpos, p.pos)

    count = zero(t)
    colobst_idx = 1
    t_to_write = 0.0

    while count < t
        # Declare these because `bt` is of un-stable type!
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
        count += increment_counter(t, tmin)
        if count > t
            count -= tmin
            break
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

    minus_d = -norm(rpos[2]-rpos[1])
    
end


function lyap_dist(p::Particle, bt::Vector{Obstacle}, t:: Float64)
    l, disp = lyapunovspectrum(p, bt, t, displacement = true)
    lyap_max = l[1]
    dist = norm(disp[2] - disp[1])

    return lyap_max, dist
end

    






##Functions that consider the observable as a pair (l,r)

function acceptance(p::Particle, s::Float64, sinphi::Float64, polygon_sides::Int, bt::Vector{Obstacle}, proposal_dist::Proposal, beta::Float64, time::Float64)
    
    l_d = x::Particle -> lyap_dist(x, bt, time)  
    lyap, r = l_d(p)
    E0 = exp(-beta*lyap)
    s1, sinphi1 = proposal_dist.f([s,sinphi])
    p1 = particle_from_coordinates([s1,sinphi1], polygon_sides, bt)
    l1, r1 = l_d(p1)
    E1 = exp(-beta*l1)
    a = min(1.0, E1/E0)
    random = rand()
    if random < a
        s = s1
        sinphi = sinphi1
        p = p1
        r = r1
        lyap = l1
    end

    return p, s, sinphi, r, lyap

end



function birkhoff_visual(t::Float64, N::Int, bt::Vector{Obstacle}, polygon_sides::Int, proposal_dist::Proposal, beta::Float64)
    
    birk_coord = zeros(N,4)
    ###initialize
    p,s, sinphi = randominside(polygon_sides, bt)
    birk_coord[1, 1:2] = [s,sinphi]
    exps, pos = lyapunovspectrum(p, bt, t; displacement = true)
    birk_coord[1, 3] = norm(pos[end] -pos[1])
    birk_coord[1,4] = exps[1]
    
    ####
    
    for i in 2:N
        p, s, sinphi, position, lambda =  acceptance(p, s, sinphi, polygon_sides, bt, proposal_dist, beta, t)
        birk_coord[i, :] = [s,sinphi, position, lambda]
    end
    
    birk_coord
end
