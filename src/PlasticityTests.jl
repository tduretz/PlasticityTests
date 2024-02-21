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
export ExtractDataTest2

include("./PlasticityModels.jl")
export MohrCoulombVermeer1990_vdev
export MohrCoulombVermeer1990_tot
export DruckerPrager_vdev

include("./StressIntegration0D_Vermeer1990.jl")
export Vermeer1990_Test1_Components
export Vermeer1990_Test1_MatVec

include("./StressIntegration0D.jl")
export Vermeer1990_StressIntegration_vdev
export Vermeer1990_StressIntegration_tot

include("./StressIntegration1D.jl")
export Main_VEP_1D
export Main_VEP_1D_tot


end # module PlasticityTests
