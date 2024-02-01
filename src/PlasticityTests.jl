module PlasticityTests

using Printf, CSV, DataFrames, Plots
using GeoParams, StaticArrays, LinearAlgebra
import LinearAlgebra
import ForwardDiff
import Statistics:mean

include("Helpers.jl")
export @minmax

include("./ReadDataFromVermeer1990.jl")
export ExtractDataCase

include("./PlasticityModels.jl")
export MohrCoulombVermeer1990_vdev
export MohrCoulombVermeer1990_tot
export DruckerPrager_vdev

include("./ElasticityModels.jl")
export ElasticModel

include("./StressIntegration0D.jl")
export Vermeer1990_StressIntegration_vdev
export Vermeer1990_StressIntegration_tot

include("./StressIntegration1D.jl")
export Main_VEP_1D


end # module PlasticityTests
