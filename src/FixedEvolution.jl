export fixed_evolution!
"""
Evolve the particle exactly until time t (slight modification of
evolve!)
"""
function fixed_evolution!(p::Particle, bt::Vector{Obstacle}, t)
    if t <= 0
        error("`evolve!()` cannot evolve backwards in time.")
    end
    rt = Float64[]
    rpos = SVector{2,Float64}[]
    rvel = SVector{2,Float64}[]
    push!(rpos, p.pos)
    push!(rvel, p.vel)
    push!(rt, 0.0)

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
            push!(rpos, p.pos + p.current_cell)
            push!(rvel, p.vel)
            push!(rt, t_to_write)
            t_to_write = 0.0
        end
    end#time loop

    tmin = t - count 
    propagate!(p, tmin)
    push!(rpos, p.pos + p.current_cell)
    push!(rvel, p.vel)
    push!(rt, t_to_write + tmin)
    

    return (rt, rpos, rvel)
end

