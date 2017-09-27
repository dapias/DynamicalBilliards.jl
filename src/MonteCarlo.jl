"""
Return the particle, its birkhoff coordinates together with the index associated with the side where the particle is located and the value of the observable
"""
function acceptance(p::Particle{T}, s::T, sinphi::T, bt::Vector{<:Obstacle{T}}, prop_dist::Proposal, cano::CanonicalDistribution, index::Int64; symmetric = true) where {T<: AbstractFloat}
    E1 = cano.f(p)
    p1, s1, sinphi1, ind = prop_dist.f(p, index)
    E2 = cano.f(p1)

    if symmetric
        if E2 - E1 >= 0.0 || cano.beta == 0.0 ##acceptance = 1.0
            p = p1
            E1 = E2
            index = ind
            s = s1
            sinphi = sinphi1
        else
            a = exp(E2 - E1)   ##recall that E = -beta*observable
            r = rand()
            if r < a   #If r is less than a the motion is accepted
                p = p1
                E1 = E2
                index = ind
            end
        end
    end

    if cano.beta !== 0.0
        observable = -E1/cano.beta
    else
        observable = E1
    end
    return p, s, sinphi, index, observable
end

"""
Performs a classical Markov Chain Monte Carlos simulation, with simmetric proposal distributions (Neighborhood and Shift)
"""

function classicalMCMC(t::T, N::Int64, bt::Vector{<:Obstacle{T}}, n::Int64, beta::T, sigma::T, tshift::T) where {T<: AbstractFloat}
    birk_coord = zeros(N,3)
    ###initialize
    p,s, sinphi, index = randominside(bt, n, sides = true)
    birk_coord[1, 1:2] = [s,sinphi]
    dist = distance(p, bt, t)
    birk_coord[1, 3] = dist    
    ####################
    ##Needed types
    obs = x::Particle -> distance(x, bt, t)
    canonical = CanonicalDistribution(beta, obs)
    shift1 = ShiftProposal(tshift, (x::Particle, i::Int64) -> shift_proposal(x, n, bt, tshift, i))
    shift2 = ShiftProposal(-tshift, (x::Particle, i::Int64) -> shift_proposal(x, n, bt, -tshift, i))
    neighborhood = NeighborhoodProposal(sigma, (x::Particle, i::Int64) -> neigborhood_proposal(x, n, bt, sigma, i) )    
    ################
    for i in 2:N
        random_proposal = rand()
        if random_proposal < 0.25
            prop_dist = shift1
        elseif 0.25 < random_proposal < 0.5
            prop_dist = shift2
        else
            prop_dist = neighborhood
        end
        p, s, sinphi, index, dist = acceptance(p, s, sinphi, bt, prop_dist, canonical, index)
         birk_coord[i, :] = [s,sinphi, dist]
    end
    
    birk_coord
end
