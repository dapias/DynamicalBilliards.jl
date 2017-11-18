__precompile__()

"""
A Julia package for dynamical billiard systems in two dimensions.

The goals of the package is to provide a flexible, easy-to-use
and intuitive framework for
fast implementation of billiard systems of arbitrary construction.
"""
module DynamicalBilliards

using StaticArrays, Distributions
import Base.show

##########################################
# Core                                   #
##########################################
include("ParticlesObstacles.jl")
include("Propagation.jl")
include("RaySplitting.jl")
include("StandardBilliards.jl")
include("LyapunovSpectrum.jl")
##For Monte Carlo simulations
include("BirkhoffMapping.jl")
include("Observables.jl")
include("MCtypes.jl")
include("MonteCarlo.jl")

##########################################
# Test Options                           #
##########################################
"""
    test_options(;print_info = false, long_tests = false)
Set if you want the long version of the tests and if you want information to be
printed during testing.
"""
function test_options(;print_info::Bool = true, long_tests::Bool = true)
    ENV["DYNAMICALBILLIARDS_PRINTTEST"] = print_info
    ENV["DYNAMICALBILLIARDS_LONGTEST"] = long_tests
end
##########################################
# Plotting Routines (loaded on demand)   #
##########################################
"""
    DynamicalBilliards.enableplotting()
Enable plotting for the package DynamicalBilliards.jl. Requires
`using PyPlot` to work properly.
"""
function enableplotting()
    dir = joinpath(dirname(@__FILE__), "plotting")
    for f in readdir(dir)
        include(joinpath(dir, f))
    end
end















end#module
