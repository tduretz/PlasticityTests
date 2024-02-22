module PlasticityTests

using Printf, CSV, DataFrames, Plots
using GeoParams, StaticArrays, LinearAlgebra
import LinearAlgebra
import ForwardDiff
import Statistics:mean

include("Helpers.jl")
export @minmax, PrincipalStress!

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

include("./StressIntegration0D_VolDev.jl")
export Vermeer1990_StressIntegration_vdev

include("./StressIntegration0D_TotalStress.jl")
export Vermeer1990_StressIntegration_tot

include("./StressIntegration1D_VolDev.jl")
export Main_VEP_1D_vdev

include("./StressIntegration1D_VolDev_GP.jl")
export Main_VEP_1D_vdev_GP

include("./StressIntegration1D_TotalStress.jl")
export Main_VEP_1D_tot

include("./StressIntegration1D_TotalStress_GP.jl")
export Main_VEP_1D_tot_GP

end # module PlasticityTests
