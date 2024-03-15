using PlasticityTests, Plots
#plotlyjs() # activate plotly backend that lets you check numbers with the mouse!

function main()

    # Case B
    σi       = (xx = -400e3, yy=-100e3, xy=0.0)

    params = (
    K    = 6.6666666667e6, # K = 3/2*Gv in Vermeer (1990)
    G    = 10e6,
    c    = 0.,
    ϕ    = 40/180*π,
    ψ    = 10/180*π,
    θt   = 25/180*π,
    ηvp  = 3e7,
    lc   = 200,
    γ̇xy  = 0.00001,
    Δt   = 20,
    nt   = 400,
    law  = :MC_Vermeer1990,
    coss = false,
    oop  = :Vermeer1990,
    pl   = true) 

    # Main_VEP_1D_vdev(σi; visu=true)
    # @time Main_VEP_1D_vdev_GP(σi; visu=false)
    # @time Main_VEP_1D_vdev_GP(σi; visu=true)

    # @time Main_VEP_1D_vdev(σi; params=params, visu=true, Ncy=20)
    @time Main_VEP_1D_vdev_coss(σi; params=params, visu=true, Ncy=20)

    # @time Main_VEP_1D_tot_GP(σi; visu=false)
    # @time Main_VEP_1D_tot(σi; visu=false)

    # @time Main_VEP_1D_vdev_coss(σi; params=params, visu=true)

end

main()