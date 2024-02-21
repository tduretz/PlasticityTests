using PlasticityTests, Plots
#plotlyjs() # activate plotly backend that lets you check numbers with the mouse!

function main()

    # Case B
    σi       = (xx = -400e3, yy=-100e3)
    # Main_VEP_1D_vdev(σi; visu=true)

    Main_VEP_1D_tot(σi; visu=true)

end

main()