export CanonicalDistribution, NeighborhoodProposal, ShiftProposal, UniformProposal, uniform_proposal, shift_proposal, neigborhood_proposal

abstract type Distribution end  #Distributions will be a type
abstract type Proposal end

"""
Type that codifies the equilibrium distribution to be sampled (the canonical one)
"""
mutable struct CanonicalDistribution{T} <: Distribution
    beta::T
    observable::Function
    f::Function
end

"""
As we know the function, we just need the parameter and the observable to construct an instance of the type
"""
CanonicalDistribution(beta, observable) =
    if beta !== 0.0
        CanonicalDistribution(beta, observable,x::Particle -> -beta*observable(x))
    else
        error("beta cannot be zero")
    end

"""
Neighborhood proposal characterized by a parameter and the proposal distribution `f`
"""
mutable struct NeighborhoodProposal{T} <: Proposal
    parameter::Union{Function, T}    #Standard_Deviation in the case of the Gaussian
    f::Function
end

"""
Shift proposal characterized by a parameter (t_shift) and the proposal distribution
"""
mutable struct ShiftProposal{T} <: Proposal
    parameter::Union{Function, T}    #t_shift
    f::Function
end

"""
The uniform proposal is constructed to pick a random point of the phase space with uniform distribution. This is codified in the function `f`
"""
mutable struct UniformProposal <: Proposal
    f::Function
end

"""
Function that is passed to the type `UniformProposal`
"""
function uniform_proposal(p::Particle{T}, n::Int, bt::Vector{<:Obstacle{T}}) where {T <: AbstractFloat}
    randominside(bt, n)
end

"""
Function that is passed to the type `NeighborhoodProposal`
"""
function neigborhood_proposal(init::InitialCondition, n::Int, bt::Vector{<:Obstacle{T}}, sigma::T) where {T <: AbstractFloat}
    
    delta_theta = rand()*(2pi)
    rprime = abs(randn()*sigma)

#    s, sinphi = coordinates_from_particle(p, n, bt, index)
    
    s_new = cos(delta_theta)*rprime + init.s
    s_new = mod(s_new, 1.0)
    sinphi_new = sin(delta_theta)*rprime + init.sinphi

    if abs(sinphi_new) > 1.0
        return init
    end

    p = particle_from_coordinates([s_new, sinphi_new], n, bt)
    side = ceil(Int, s_new*n)

    return InitialCondition(p, s_new, sinphi_new, side)
end

"""
Function that evolves a particle until it reaches a periodic wall. It returns the index associated with the periodic image of the side that the particle reaches (where the motion would continue)
"""
function evolution_to_periodicwall!(p::Particle{T}, bt::Vector{<:Obstacle{T}}, n::Int) where {T<: AbstractFloat}
    i = 1
    while true
        tmin::T, i::Int = next_collision(p, bt)
        tmin = DynamicalBilliards.relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        if typeof(bt[i]) <: PeriodicWall
            break
        end
    end
    if n == 6
        index = i >= 4 ? mod(i + div(n,2), n) : i + div(n,2)
    elseif n==4
        index = i >= 3 ? mod(i + div(n,2), n) : i + div(n,2)
    end
        
end

"""
Function that is passed to the type `ShiftProposal`
"""
function shift_proposal(particle::Particle{T}, sides::Int, bt::Vector{<:Obstacle{T}}, tshift::T, index::Int) where {T <: AbstractFloat}
    p = copy(particle)
    if tshift < 0.0  ##Invert the direction of the velocity
        #and put the particle on its periodic image
        p.vel = -p.vel
        p.pos += bt[index].normal 
    end
    
    t = abs(tshift)
    count = zero(T)
    
    while count < t
        tmin::T, i::Int = next_collision(p, bt)
        # set counter
        if count +  DynamicalBilliards.increment_counter(t, tmin) > t
            break
        end
        #####
        tmin = DynamicalBilliards.relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        count += DynamicalBilliards.increment_counter(t, tmin)
    end#time loop
    ##Here the evolution ends and then we propagate until the next collision
    new_index = evolution_to_periodicwall!(p, bt, sides)
    new_particle = Particle(p.pos, p.vel, SVector{2,T}([0.,0.]))

    if tshift < 0.0
        new_particle.vel = -new_particle.vel
        new_particle.pos += bt[new_index].normal
        if sides == 6
            new_index = new_index >= 4 ? mod(new_index + div(sides,2), sides) : new_index + div(sides,2)
        elseif sides == 4
            new_index = new_index >= 3 ? mod(new_index + div(sides,2), sides) : new_index + div(sides,2)
        end
    end
    s, sinphi = coordinates_from_particle(new_particle, sides, bt, new_index)
    
    return InitialCondition(new_particle, s, sinphi, new_index)
end


