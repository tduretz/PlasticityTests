using PlasticityTests, Plots, LinearAlgebra

function main()

    # MC params
    ϕ  = 40.0*π/180.
    ψ  = 10.0*π/180.
    G  = 10e6
    c  = 0.

    𝐃ᵉ  = [2G 0 0; 0 2G 0; 0 0 G]

    Δt  = 50
    nt  = 400
    γ̇xy = 0.00001

    θ  = π/2
    θₐ = π/4 + 0.25*(ϕ + ψ)

    σh = -100e3
    σv = (1 + sin(ϕ)) / (1 - sin(ϕ)) * σh

    # to Cartesian
    σxx_o = 1/2*(σh + σv) +  1/2*(σh - σv)*cos(2*θₐ)
    σyy_o = 1/2*(σh + σv) -  1/2*(σh - σv)*cos(2*θₐ)
    σxy_o = 1/2*(σh - σv)*sin(2*θₐ)

    σxx_i = 1/2*(σh + σv) +  1/2*(σh - σv)*cos(2*θₐ)
    σyy_i = 1/2*(σh + σv) -  1/2*(σh - σv)*cos(2*θₐ)
    σxy_i = 1/2*(σh - σv)*sin(2*θₐ)
    σ_o   = [σxx_o; σyy_o; σxy_o]
    σ_i   = [σxx_i; σyy_i; σxy_i]

    @show θₐ
    @show 0.5*atan(2*σxy_i/(σxx_i-σyy_i))

    Δγxy = γ̇xy*Δt

    load = zeros(nt)

    for it=1:1#nt

        # Stress outside
        σ_o  = [σxx_o; σyy_o; σxy_o]
        τ_o  = sqrt(0.25*(σ_o[1]-σ_o[2])^2 + σ_o[3]^2)
        f    = τ_o + 0.5*(σ_o[1]+σ_o[2])*sin(ϕ) 
        θ_o  = 0.5*acos((σ_o[1]-σ_o[2])/2/τ_o)
        θ    = θ_o

        # Stress integration
        𝐚    = G.*[ cos(2θ) + sin(ψ)
                -cos(2θ) + sin(ψ) 
                    sin(2θ)         ] 
        𝐛    = G.*[ cos(2θ) + sin(ϕ)
                -cos(2θ) + sin(ϕ) 
                    sin(2θ)         ] 
        d    = G*(1 + sin(ψ)*sin(ϕ))
        𝐌    = 𝐃ᵉ - (f>=0)*1/d*𝐚*𝐛'
        
        α    = 1/2*(1 + cos(2*θₐ) - (1 - cos(2*θₐ)*𝐌[2,1]/2/G))
        β    = 1/2*(  - sin(2*θ) - (1 - cos(2*θₐ)*𝐌[3,1]/2/G))
        Δσv  = G^2*(cos(2*θ) + sin(ψ)) * (cos(2*θ) + sin(ϕ)) / ( (1 + sin(ϕ)*sin(ψ)) * (β*𝐌[2,2] + α*𝐌[3,2])) * Δγxy

        σv += Δσv 

        # τ_o    = sqrt(0.25*(σ_o[1]-σ_o[2])^2 + σ_o[3]^2)
        # θ_o    = 0.5*acos((σ_o[1]-σ_o[2])/2/τ_o)
        σxx_o  = 1/2*(σh + σv) +  1/2*(σh - σv)*cos(2*θₐ)
        σyy_o  = 1/2*(σh + σv) -  1/2*(σh - σv)*cos(2*θₐ)
        σxy_o  = 1/2*(σh - σv)*sin(2*θₐ)

        σxx_i  = σxx_o
        σyy_i  = σyy_o
        σxy_i  = σxy_o
        σ_i    = [σxx_i; σyy_o; σxy_o]
        τ_i    = sqrt(0.25*(σ_i[1]-σ_i[2])^2 + σ_i[3]^2)
        f      = τ_i + 0.5*(σ_i[1]+σ_i[2])*sin(ϕ)
        @show fc     = f
        if fc>0
            for iter=1:500
                σ_i   .= [σxx_i; σyy_o; σxy_o]
                τ_i    = sqrt(0.25*(σ_i[1]-σ_i[2])^2 + σ_i[3]^2)
                fc     = τ_i + 0.5*(σ_i[1]+σ_i[2])*sin(ϕ)
                σxx_i -= fc
                if abs(fc)<1e-8 break end
            end
        end
        @show (f,fc)
        σ_i  = [σxx_i; σyy_o; σxy_o]
        τ_i  = sqrt(0.25*(σ_i[1]-σ_i[2])^2 + σ_i[3]^2)
        θ_i  = 0.5*acos((σ_i[1]-σ_i[2])/2/τ_i)

        # MC out
        σzz_o = 1/2*(σxx_o + σyy_o)
        σ     = [σxx_o σxy_o 0; σxy_o σyy_o 0; 0 0 σzz_o]
        sp    = eigvals(σ) 
        σ1_o  = sp[1]  
        σ3_o  = sp[3]  

        # MC in
        σzz_i = 1/2*(σxx_i + σyy_i)
        σ     = [σxx_i σxy_i 0; σxy_i σyy_i 0; 0 0 σzz_i]
        sp    = eigvals(σ) 
        σ1_i  = sp[1]  
        σ3_i  = sp[3]

        θmc   = LinRange(-π, 0, 100)
        σMC   = LinRange(-500, 0, 100 ) .*1e3
        τMC   = -σMC.*tan(ϕ) .+ c
        yield = (x = σMC./1e3, y = τMC./1e3)
        τA    = (σ1_o - σ3_o)/2
        PA    = (σ1_o + σ3_o)/2
        τB    = (σ1_i - σ3_i)/2
        PB    = (σ1_i + σ3_i)/2
        MC_A  = (x = (τA.*cos.(θmc) .+ PA)./1e3, y = (τA.*sin.(θmc))./1e3) 
        MC_B  = (x = (τB.*cos.(θmc) .+ PB)./1e3, y = (τB.*sin.(θmc))./1e3)

        load[it] = σv/σh

        # if mod(it,50)==0

        p1 = plot(title="Test 1 Mohr circles", ylabel="τ", xlabel="σₙ", size=(300,300), aspect_ratio=1)
        p1 = plot!( MC_A... , color=:blue, label="Case A" )
        p1 = plot!( MC_B...,  color=:green, label="Case B"  )
        p1 = plot!( yield..., color=:red, label="Yield"  )

        p2 = plot((1:it).*Δt.*γ̇xy*100, load[1:it])


        display(plot(p1, p2))
        sleep(0.1)
        # end
    end
end

main()