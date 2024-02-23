using PlasticityTests, Plots
#plotlyjs() # activate plotly backend that lets you check numbers with the mouse!

function main()

    # Case B
    σi       = (xx = -400e3, yy=-100e3, xy=0.0)
    # Main_VEP_1D_vdev(σi; visu=true)
    # @time Main_VEP_1D_vdev_GP(σi; visu=true)
    
    # @time Main_VEP_1D_tot_GP(σi; visu=true)
    # @time Main_VEP_1D_tot(σi; visu=true)

    @time Main_VEP_1D_vdev_coss(σi; visu=true)

end

main()