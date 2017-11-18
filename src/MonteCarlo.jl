export  chaoticMCMC, t_star, lambda_tstar, sigma_local, symmetricMCMC, uniformMCMC, directsampling, ratio_proposals, acceptance, newMCMC


"""
Equation 54, Leitao, Lopes, Altmann
"""
function t_star(p::Particle{T}, bt::Vector{<:Obstacle{T}}, beta::T, to::T, lambda_L::T; a = 0.5) where {T <: AbstractFloat}
    lambda_to = lyapunovmaximum(p, bt, to) 
    tstar = to - abs((a-1.)/beta * 1/(lambda_L - lambda_to))
    if tstar < 0.0
        return T(0.0)
    end
    tstar
end

function t_star(p::Particle{T}, bt::Vector{<:Obstacle{T}}, beta::T, to::T, D::T, observable::Function; a = 0.5) where {T <: AbstractFloat}
    y = observable(p)
    tstar = to - abs((a-1.)/(2*D*beta) * 1/(1 - y^2.))
#    A = abs((1-a)/(2*D*beta))
#    num = -A + to + A*y^2 - 2*to*y^2 + 3*to*y^4 - 2*sqrt(-A^2*y^4 + A*to*y^4 - A*to*y^6 + to^2*y^8)
    #    den = 1 - 2*y^2 + 5*y^4
#    num = -2A + 2to + 2A*y^2 - 4to*y^2 + 3*to*y^4 - sqrt(-4*A^2*y^4 + 4*A*to*y^4 - 4*A*to*y^6 + to^2*y^8)
#    den = 2*(1 - 2y^2 + 2y^4)
#    tstar = num/den
    if tstar < 0.0
        return T(0.0)
    end
    tstar
end
"""
Lyapunov exponent at time tstar
"""
function lambda_tstar(p::Particle{T}, bt::Vector{<:Obstacle{T}}, tstar::T)  where {T <: AbstractFloat}
    if tstar > 0.0
        return lyapunovmaximum(p, bt, tstar) 
    else
        return 0.0
    end
end
"""
Equation 43, Leitao, Lopes, Altmann
"""
function sigma_local(p::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, lambda_L::T; Delta = T(1.0))  where {T <: AbstractFloat}
    tstar =  t_star(p, bt, beta, t, lambda_L)
    l_tstar = lambda_tstar(p, bt, tstar)
    sigma = Delta*exp(-l_tstar*tstar)
end

function sigma_local(p::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, D::T, ysq::T; Delta = T(1.0))  where {T <: AbstractFloat}
    tstar =  t_star(p, bt, beta, t, D, ysq)
    l_tstar = T(2.10)
    sigma = Delta*exp(-l_tstar*tstar)
end

function t_star(p::Particle{T}, bt::Vector{<:Obstacle{T}}, beta::T, to::T, D::T, ysq::T; a = 0.5) where {T <: AbstractFloat}
    tstar = to - abs((a-1.)/(2*D*beta) * 1/(1 - ysq))
    if tstar < 0.0
        return T(0.0)
    end
    tstar
end

function sigma_local(p::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, D::T, observable::Function; Delta = T(1.0))  where {T <: AbstractFloat}
    tstar =  t_star(p, bt, beta, t, D, observable)
#    l_tstar = lambda_tstar(p, bt, tstar)
    l_tstar = T(2.10)
    sigma = Delta*exp(-l_tstar*tstar)
end



"""
Compute the distance between two points on the set of initial conditions
"""
function distance(i1::InitialCondition{T}, i2::InitialCondition{T}) where {T<:AbstractFloat}
    norm([i1.s - i2.s, i1.sinphi - i2.sinphi])
end
"""
Compute the factor g(x' --> x)/g(x --> x') 
"""
function ratio_proposals(prop::NeighborhoodProposal, init1::InitialCondition{T}, init2::InitialCondition{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, lambda_L::T) where {T<:AbstractFloat}
    d  = distance(init1, init2)
    sigma1 = prop.parameter  
    sigma2 = sigma_local(init2.particle, bt, t, beta, lambda_L)
    expo = -d^2/2*(1./sigma2^2. - 1./sigma1^2.)
    return exp(expo)*sigma1/sigma2
end

function ratio_proposals(prop::NeighborhoodProposal, init1::InitialCondition{T}, init2::InitialCondition{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, D::T, obs::Function; delta = T(1.0)) where {T<:AbstractFloat}
    d  = distance(init1, init2)
    sigma1 = prop.parameter  
    sigma2 = sigma_local(init2.particle, bt, t, beta, D, obs, Delta = delta)
    expo = -d^2/2*(1./sigma2^2. - 1./sigma1^2.)
    return exp(expo)*sigma1/sigma2
end

function ratio_proposals(prop::NeighborhoodProposal, init1::InitialCondition{T}, init2::InitialCondition{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, D::T, ysq::T; delta = T(1.0)) where {T<:AbstractFloat}
    d  = distance(init1, init2)
    sigma1 = prop.parameter  
    sigma2 = sigma_local(init2.particle, bt, t, beta, D, ysq)
    expo = -d^2/2*(1./sigma2^2. - 1./sigma1^2.)
    return exp(expo)*sigma1/sigma2
end


function ratio_proposals(prop::ShiftProposal, init1::InitialCondition{T}, init2::InitialCondition{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, lambda_L::T) where {T<:AbstractFloat}
    return  T(1.0)
end

function ratio_proposals(prop::ShiftProposal, init1::InitialCondition{T}, init2::InitialCondition{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, D::T, obs::Function) where {T<:AbstractFloat}
    return  T(1.0)
end

function ratio_proposals(prop::ShiftProposal, init1::InitialCondition{T}, init2::InitialCondition{T}, bt::Vector{<:Obstacle{T}}, t::T, beta::T, D::T, ysq::T) where {T<:AbstractFloat}
    return  T(1.0)
end

"""
Return the particle, its birkhoff coordinates together with the index associated with the side where the particle is located and the value of the observable
"""
function acceptance(init1::InitialCondition, bt::Vector{<:Obstacle{T}}, prop_dist::Proposal, cano::CanonicalDistribution; symmetric = true,  lambda_mean= 1.0, time = 1.0) where {T<: AbstractFloat}
    E1 = cano.f(init1.particle)
    init2 = prop_dist.f(init1)
    E2 = cano.f(init2.particle)
    init = init1

    if symmetric
        acceptance = exp(E2 - E1)   ##recall that E = -beta*observable
        r = rand()
        if r < acceptance  #If r is less than a the motion is accepted
            init = init2
            E1 = E2
        end

    else        
        ratio = ratio_proposals(prop_dist, init1, init2, bt, time, cano.beta, lambda_mean)
        acceptance = ratio*exp(E2 - E1)
        r = rand()
        if r < acceptance   #If r is less than a the motion is accepted
            init = init2
            E1 = E2
        end
    end

    observable = -E1/cano.beta

    return init, observable
end

"""
Performs a classical Markov Chain Monte Carlos simulation, with simmetric proposal distributions (Neighborhood and Shift)
"""
function symmetricMCMC(t::T, N::Int, bt::Vector{<:Obstacle{T}}, n::Int, beta::T, sigma::T, tshift::T; D = T(0.17)) where {T<: AbstractFloat}
    birk_coord = zeros(N,3)
    ###initialize
    obs = x::Particle -> distance(x, bt, t)/sqrt(2*D*t)
    init  = randominside(bt, n)
    birk_coord[1, 1:2] = [init.s, init.sinphi]
    dist = obs(init.particle)
    birk_coord[1, 3] = dist    
    ####################
    ##Needed types
    canonical = CanonicalDistribution(beta, obs)
    shift1 = ShiftProposal(tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, tshift, x.index))
    shift2 = ShiftProposal(-tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, -tshift, x.index))
    neighborhood = NeighborhoodProposal(sigma, x::InitialCondition -> neigborhood_proposal(x, n, bt, sigma) )    
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
        init, dist = acceptance(init, bt, prop_dist, canonical)
        birk_coord[i, :] = [init.s,init.sinphi, dist]
    end
    
    birk_coord
end
"""
Performs a simulation based on a uniform sampling of the initial set using Birkhoff coordinates
"""
function directsampling(t::T, N::Int, bt::Vector{<:Obstacle{T}}, n::Int; lyapunov= false, observable=distance) where {T<: AbstractFloat}
    if lyapunov
        birk_coord = zeros(T,N,4)
    else
        birk_coord = zeros(T, N,3)
    end       
    ###initialize
    init  = randominside(bt, n)
    birk_coord[1, 1:2] = [init.s, init.sinphi]
    dist = observable(init.particle, bt, t)
    birk_coord[1, 3] = dist
    if lyapunov
        birk_coord[1,4] = lyapunovmaximum(init.particle, bt, t)
    end
    ####################
    ##Needed types
    obs = x::Particle -> observable(x, bt, t)
    lyap = x::Particle -> lyapunovmaximum(x, bt,t)
    uniform  = UniformProposal(x::Particle -> uniform_proposal(x, n, bt))
    ################
    for i in 2:N
        init = uniform.f(init.particle)
        dist = obs(init.particle)
        if lyapunov
            l = lyap(init.particle)
            birk_coord[i, :] = [init.s, init.sinphi, dist, l]
        else
            birk_coord[i, :] = [init.s, init.sinphi, dist]
        end
    end
    
    birk_coord
end
"""
Function that implements the algorithm set up in the article *Importance Sampling of Rare Events in Chaotic Systems*
"""
function chaoticMCMC(t::T, N::Int, bt::Vector{<:Obstacle{T}}, n::Int, beta::T, tshift::T, lambda_L::T) where {T<: AbstractFloat}
    birk_coord = zeros(N,3)
    ###initialize
    init  = randominside(bt, n)
    birk_coord[1, 1:2] = [init.s, init.sinphi]
    dist = distance(init.particle, bt, t)
    birk_coord[1, 3] = dist    
    ####################
    ##Needed types
    obs = x::Particle -> distance(x, bt, t)
    canonical = CanonicalDistribution(beta, obs)
    shift1 = ShiftProposal(tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, tshift, x.index))
    shift2 = ShiftProposal(-tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, -tshift, x.index))
    ################
    for i in 2:N
        random_proposal = rand()
        if random_proposal < 0.25
            prop_dist = shift1
        elseif 0.25 < random_proposal < 0.5
            prop_dist = shift2
        else
            sigma = sigma_local(init.particle, bt, t, beta, lambda_L)
            prop_dist = NeighborhoodProposal(sigma, x::InitialCondition -> neigborhood_proposal(x, n, bt, sigma) )
        end
        init, dist = acceptance(init, bt, prop_dist, canonical, symmetric = false, lambda_mean = lambda_L, time = t)
         birk_coord[i, :] = [init.s, init.sinphi, dist]
    end
    
    birk_coord
end
"""
Classical Importance Sampling using a uniform proposal
"""
function uniformMCMC(t::T, N::Int, bt::Vector{<:Obstacle{T}}, n::Int, beta::T) where {T<: AbstractFloat}
    birk_coord = zeros(N,3)
    ###initialize
    init  = randominside(bt, n)
    birk_coord[1, 1:2] = [init.s, init.sinphi]
    dist = distance(init.particle, bt, t)
    birk_coord[1, 3] = dist    
    ####################
    ##Needed types
    obs = x::Particle -> distance(x, bt, t)
    ################
    dist = obs(init.particle)
    for i in 2:N
        init2 = randominside(bt, n)
        dist2 = obs(init2.particle)
        acceptance = exp(-beta*(dist2 - dist))
        r = rand()
        if r < acceptance  #If r is less than a the motion is accepted
            init = init2
            dist = dist2
        end
        birk_coord[i, :] = [init.s,init.sinphi, dist]
    end
    
    birk_coord
end




function newMCMC(t::T, N::Int, bt::Vector{<:Obstacle{T}}, n::Int, beta::T, tshift::T, D::T) where {T<: AbstractFloat}
    birk_coord = zeros(N,3)
    ###initialize
    obs = x::Particle -> distance(x, bt, t)/sqrt(2*D*t)
    init  = randominside(bt, n)
    birk_coord[1, 1:2] = [init.s, init.sinphi]
    dist = obs(init.particle)  ##Actually this is y
    birk_coord[1, 3] = dist    
    ####################
    ##Needed types
    canonical = CanonicalDistribution(beta, obs)
    shift1 = ShiftProposal(tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, tshift, x.index))
    shift2 = ShiftProposal(-tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, -tshift, x.index))
    ################
    for i in 2:N
        random_proposal = rand()
        if random_proposal < 0.25
            prop_dist = shift1
        elseif 0.25 < random_proposal < 0.5
            prop_dist = shift2
        else
            sigma = sigma_local(init.particle, bt, t, beta, D, obs)
            prop_dist = NeighborhoodProposal(sigma, x::InitialCondition -> neigborhood_proposal(x, n, bt, sigma) )
        end
        init, dist = acceptance(init, bt, prop_dist, canonical, D, obs, t)
        birk_coord[i, :] = [init.s, init.sinphi, dist]
    end
    
    birk_coord
end

function newMCMC(t::T, N::Int, bt::Vector{<:Obstacle{T}}, n::Int, beta::T, tshift::T, D::T, obs::Function) where {T<: AbstractFloat}
    birk_coord = zeros(N,3)
    ###initialize
    init  = randominside(bt, n)
    birk_coord[1, 1:2] = [init.s, init.sinphi]
    dist = obs(init.particle)  ##Actually this is y
    birk_coord[1, 3] = dist    
    ####################
    ##Needed types
    canonical = CanonicalDistribution(beta, obs)
    shift1 = ShiftProposal(tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, tshift, x.index))
    shift2 = ShiftProposal(-tshift, x::InitialCondition -> shift_proposal(x.particle, n, bt, -tshift, x.index))
    ################
    for i in 2:N
        random_proposal = rand()
        ysq = obs(init.particle)
        if random_proposal < 0.25
            prop_dist = shift1
        elseif 0.25 < random_proposal < 0.5
            prop_dist = shift2
        else
            sigma = sigma_local(init.particle, bt, t, beta, D, ysq)
            prop_dist = NeighborhoodProposal(sigma, x::InitialCondition -> neigborhood_proposal(x, n, bt, sigma) )
        end
        init, dist = acceptance(init, bt, prop_dist, canonical, D, ysq, t)
        birk_coord[i, :] = [init.s, init.sinphi, dist]
    end
    
    birk_coord
end


function acceptance(init1::InitialCondition, bt::Vector{<:Obstacle{T}}, prop_dist::Proposal, cano::CanonicalDistribution, D::T, obs::Function, time::T) where {T<: AbstractFloat}
    E1 = cano.f(init1.particle)
    init2 = prop_dist.f(init1)
    E2 = cano.f(init2.particle)
    init = init1
    ratio = ratio_proposals(prop_dist, init1, init2, bt, time, cano.beta, D, obs)
    acceptance = ratio*exp(E2 - E1)
    r = rand()
    if r < acceptance   #If r is less than a the motion is accepted
        init = init2
        E1 = E2
    end
    
    observable = -E1/cano.beta

    return init, observable
end

function acceptance(init1::InitialCondition, bt::Vector{<:Obstacle{T}}, prop_dist::Proposal, cano::CanonicalDistribution, D::T, ysq::T, time::T) where {T<: AbstractFloat}
    E1 = cano.f(init1.particle)
    init2 = prop_dist.f(init1)
    E2 = cano.f(init2.particle)
    init = init1
    ratio = ratio_proposals(prop_dist, init1, init2, bt, time, cano.beta, D, ysq)
    acceptance = ratio*exp(E2 - E1)
    r = rand()
    if r < acceptance   #If r is less than a the motion is accepted
        init = init2
        E1 = E2
    end
    
    observable = -E1/cano.beta

    return init, observable
end


