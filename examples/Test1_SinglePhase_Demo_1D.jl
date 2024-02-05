using PlasticityTests, Plots

function main()

    # Case B
    σi       = (xx = -400e3, yy=-100e3)
    Main_VEP_1D(σi; visu=true)

end

main()