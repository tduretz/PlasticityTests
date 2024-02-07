using Plots, LinearAlgebra, CSV

# Vermeer (1990)
# cohesion is 0!
# ν is 0 !
# Κ is 2/3G
# σzzi = -p
# Vermeer integration formulae: should include elasticity in the first steps!
# Total OOP strain is 1/2*(Exx+Eyy)
# Total OOP plastic strain is 0
# Deviatoric OOP plastic strain is -1/2*(exxp+eyyp)

function Vermeer1990_FirstTest_Components(σxxi, σyyi, params)

    # Material properties
    G   = params.G
    ϕ   = params.ϕ
    ψ   = params.ψ
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    # Initial condition
    σxx = σxxi 
    σyy = σyyi
    σxy = 0. 
    σm  = [σxx σxy; σxy σyy]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
    
    # Storage
    app_fric    = zeros(nt)
    γxy         = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σxy/σyy

    for it=2:nt

        # σ3 angle
        σm   = [σxx σxy; σxy σyy]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
        # θ    =  0.5*asin(σxy/τ)
        θ    = 0.5*acos((σxx-σyy)/2/τ)

        # Stress integration
        Gs   = G*(2*sin(2θ)^2 + (cos(2θ) + sin(ψ))*(cos(2θ) + sin(ϕ)))^(-1)
        σ̇xx  = -2Gs*sin(2θ)*(cos(2θ) + sin(ψ))*γ̇xy
        σ̇yy  = 0.0
        σ̇xy  = Gs*(cos(2θ) + sin(ψ))*(cos(2θ) + sin(ϕ))*γ̇xy
        σxx += σ̇xx*Δt 
        σyy += σ̇yy*Δt
        σxy += σ̇xy*Δt

        # Storage
        app_fric[it] = -σxy/σyy
        γxy[it]      = γxy[it-1] + γ̇xy*Δt
        θv[it]       = θ
    end
    return γxy, θv, app_fric
end

function Vermeer1990_FirstTest_Matrix(σxxi, σyyi, params)

    # Material properties
    G   = params.G
    ϕ   = params.ϕ
    ψ   = params.ψ
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    𝐃ᵉ  = [2G 0 0; 0 2G 0; 0 0 G]
    
    # Initial condition
    σ   = [σxxi; σyyi; 0.]
    σm  = [σ[1] σ[3]; σ[3] σ[2]]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
 
    # Storage
    app_fric    = zeros(nt)
    γxy         = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σ[3]/σ[2]

    for it=2:nt

        # σ3 angle
        σm   = [σ[1] σ[3]; σ[3] σ[2]]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σ[1]-σ[2])^2 + σ[3]^2)
        θ    = 0.5*acos((σ[1]-σ[2])/2/τ)
        
        f      = τ + 0.5*(σ[1]+σ[2])*sin(ϕ) 

        # Stress integration
        𝐚    = G.*[ cos(2θ) + sin(ψ)
                   -cos(2θ) + sin(ψ) 
                    sin(2θ)         ] 
        𝐛    = G.*[ cos(2θ) + sin(ϕ)
                -cos(2θ) + sin(ϕ) 
                 sin(2θ)         ] 
        d    = G*(1 + sin(ψ)*sin(ϕ))
        𝐌    = 𝐃ᵉ - (f>=0)*1/d*𝐚*𝐛'
        ε̇    = [0.;   -(𝐌[2,3]/𝐌[2,2])*γ̇xy;   γ̇xy]
        σ  .+= 𝐌*ε̇*Δt

        τ    = sqrt(0.25*(σ[1]-σ[2])^2 + σ[3]^2)
        θ    = 0.5*acos((σ[1]-σ[2])/2/τ)
        
        # Storage
        app_fric[it] = -σ[3]/σ[2]
        γxy[it]      = γxy[it-1] + γ̇xy*Δt
        θv[it]       = θ
    end
    return γxy, θv, app_fric
end

function Vermeer1990_StressIntegration1(σxxi, σyyi, params)

    @info "Stress integration: total stress"

    # Material properties
    G   = params.G
    ϕ   = params.ϕ
    ψ   = params.ψ
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    # Initial condition
    σxx = σxxi 
    σyy = σyyi
    σxy = 0. 
    σzz = 0.5*(σxxi+σyyi)
    σm  = [σxx σxy; σxy σyy]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
    
    # Storage
    app_fric    = zeros(nt)
    γxy         = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σxy/σyy

    ε̇xx = 0.
    ε̇yy = 0. # to be found
 
    for it=2:nt

        σxx0 = σxx
        σyy0 = σyy
        σzz0 = σzz
        σxy0 = σxy
        λ = 0.

        # σ3 angle
        σm   = [σxx σxy; σxy σyy]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
        # θ    =  0.5*asin(-σxy/τ) + π/2
        # θ    = 0.5*acos((σxx-σyy)/2/τ)

        # Trial stress integration
        ε̇yy = (σyy0 - σyyi) / (2*G*Δt)
        σxx = σxx0 + 2*G*Δt*ε̇xx
        σyy = σyy0 + 2*G*Δt*ε̇yy # this will not change
        σxy = σxy0 +   G*Δt*γ̇xy
        ε̇zz = 0. + (σzz0 - σzz) / (2*G*Δt)
        ε̇xxd = ε̇xx - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        ε̇yyd = ε̇yy - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        ε̇zzd = ε̇zz - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        σzz =  (σxx + σyy)/2
        P   = -(σxx + σyy)/2

        # Yield 
        τ      = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
        f      = τ + 0.5*(σxx+σyy)*sin(ϕ) 
        fc     = f
        if f>0 && params.pl
            dQdσxx = 0.25*(σxx-σyy)/τ + 0.5*sin(ψ)
            dQdσyy =-0.25*(σxx-σyy)/τ + 0.5*sin(ψ)
            dQdσxy = σxy/τ
            λ = 0.0
            for iter=1:100
                ε̇xxp = λ*dQdσxx
                ε̇yyp = λ*dQdσyy
                # ε̇zzp = 1/2*(ε̇xxp + ε̇yyp)
                ε̇zzp = 0.0
                γ̇xyp = λ*dQdσxy
                ε̇yy  = ε̇yyp + (σyy0 - σyyi) / (2*G*Δt)
                ε̇zz  = 1/2*(ε̇xx + ε̇yy)
                σxx  = σxx0 + 2*G*Δt*(ε̇xx - ε̇xxp)
                σyy  = σyy0 + 2*G*Δt*(ε̇yy - ε̇yyp) # this will not change
                σxy  = σxy0 +   G*Δt*(γ̇xy - γ̇xyp)
                σzz  = σzz0 + 2*G*Δt*(ε̇zz - ε̇zzp)
                ε̇xxd = ε̇xx - 1/3*(ε̇xx + ε̇yy + ε̇zz)
                ε̇yyd = ε̇yy - 1/3*(ε̇xx + ε̇yy + ε̇zz)
                ε̇zzd = ε̇zz - 1/3*(ε̇xx + ε̇yy + ε̇zz)
                P    = -1/3*(σxx + σyy + σzz)
                τ    = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
                fc   = τ + 0.5*(σxx+σyy)*sin(ϕ)
                λ   += fc / G / 10
                if abs(f)<1e-8 break end
            end
        end

        if it==30 
            @printf("%d %f %2.2e %f %2.2e %f %f %2.2e %2.2e %2.2e\n", it, f, λ, P, ε̇yy, σyy, σzz+P, ε̇zz, ε̇zzd, ε̇xx + ε̇yy + ε̇zz)
            # @show  (σxx + σyy)/2, σzz
            # @show -(σxx + σyy)/2, P
            τxx = σxx+P
            τyy = σyy+P
            τzz = σzz+P
            @printf("%1.2f %1.2f %1.2f %1.2f %1.2f\n", σxx, σyy, σzz, σxy, -1/3*(σxx+σyy+σzz))
            @show ε̇xx
            @show ε̇yy
            @show ε̇zz
            @show ε̇zzd
            @show τxx
            @show τyy
            @show τzz
            @show σzz
            @show -1/3*(τxx+τyy+τzz-3P), P
            @show τxx+τyy+τzz
            @show  ε̇zzd + ε̇xxd + ε̇yyd
        end
        
        # Storage
        app_fric[it] = -σxy/σyy
        γxy[it]      = γxy[it-1] + γ̇xy*Δt
        θv[it]       = θ
    end
    return γxy, θv, app_fric
end

function Vermeer1990_StressIntegration2(σxxi, σyyi, params)

    @info "Stress integration: volumetric-deviatoric split"

    # Material properties
    G   = params.G
    ϕ   = params.ϕ
    ψ   = params.ψ
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    # Initial condition
    P    = -1/2*(σxxi+σyyi)  # σzzi = -p
    τxxi = σxxi + P
    τyyi = σyyi + P
    τxx = τxxi 
    τyy = τyyi
    τzz = 0.
    τxy = 0. 

    σxx = σxxi 
    σyy = σyyi
    σxy = 0. 
    σzz = 0.5*(σxxi+σyyi)
    
    σm  = [τxx-P τxy; τxy τyy-P]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
    
    # Storage
    app_fric    = zeros(nt)
    γxy         = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -τxy/(τyy-P)
 
    for it=2:nt

        τxx0 = τxx
        τyy0 = τyy
        τzz0 = τzz
        τxy0 = τxy
        P0   = P

        σxx0  = σxx
        σyy0  = σyy
        σxy0  = σxy
        σzz0  = σzz

        λ = 0.

        # σ3 angle
        σm  = [τxx-P τxy; τxy τyy-P]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
     
        # Trial stress integration
        ε̇xxp, ε̇yyp, ε̇zzp, γ̇xyp, ∇vp = 0., 0., 0., 0., 0.
        ε̇xx = 0.
        ε̇yy = 0. # to be found
        ε̇yy = (0.333333333333333 * G .* Δt .* (-ε̇xxp + 2.0 * ε̇yyp - ε̇zzp + ∇vp) + 0.5 * P0 + 0.5 * σyyi - 0.5 * τyy0) ./ (G .* Δt)
        
        # ε̇yy = 0.166666666666667 * (G .* Δt .* (-3.0 * ε̇xxp + 3.0 * ε̇yyp + 2.0 * ∇vp) + 3.0 * P0 + 3.0 * σyyi - 3.0 * τyy0) ./ (G .* Δt)
        ε̇zz = 1/2*(ε̇xx + ε̇yy)
        ∇v  = ε̇xx + ε̇yy + ε̇zz
        τxx = τxx0 + 4/3*G*Δt*ε̇xx - 2/3*G*Δt*ε̇yy - 2/3*G*Δt*ε̇zz
        τyy = τyy0 + 4/3*G*Δt*ε̇yy - 2/3*G*Δt*ε̇xx - 2/3*G*Δt*ε̇zz
        τzz = τzz0 + 4/3*G*Δt*ε̇zz - 2/3*G*Δt*ε̇xx - 2/3*G*Δt*ε̇yy
        τxy = τxy0 + G*Δt*γ̇xy
        P   = P0   - 2/3*G*Δt*∇v

        ε̇xxd = ε̇xx - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        ε̇yyd = ε̇yy - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        ε̇zzd = ε̇zz - 1/3*(ε̇xx + ε̇yy + ε̇zz)

        σxx  = σxx0 + 2*G*Δt*(ε̇xx - ε̇xxp)
        σyy  = σyy0 + 2*G*Δt*(ε̇yy - ε̇yyp) # this will not change
        σxy  = σxy0 +   G*Δt*(γ̇xy - γ̇xyp)
        σzz  = σzz0 + 2*G*Δt*(ε̇zz - ε̇zzp)

        # Yield 
        τ    = sqrt(0.25*(τxx-τyy)^2 + τxy^2)
        f    = τ + 0.5*(τxx+τyy-2P)*sin(ϕ) 
        fc   = f
        if f>0 && params.pl
            σxx = τxx - P
            σyy = τyy - P
            dQdσxx = 0.25*(σxx-σyy)/τ + 0.5*sin(ψ)
            dQdσyy =-0.25*(σxx-σyy)/τ + 0.5*sin(ψ)
            dQdσzz = 0.
            dQdτxx =  0.25*(τxx-τyy)/τ + 0.5*sin(ψ)
            dQdτyy = -0.25*(τxx-τyy)/τ + 0.5*sin(ψ)
            dQdτzz = 0
            dQdτxy =  τxy/τ
            dQdP   = -sin(ψ)
            dQdτxx = dQdσxx - 0/3*(dQdσxx+dQdσyy)
            dQdτyy = dQdσyy - 0/3*(dQdσxx+dQdσyy)
            λ      = 0.0
            for iter=1:100
                γ̇xyp =  λ*dQdτxy
                ε̇xxp =  λ*dQdτxx 
                ε̇yyp =  λ*dQdτyy

                # with this one: ezzp'=0, ezzp!=0, tzz=0, ezz'=0
                # ε̇zzp =  1/2*(ε̇xxp + ε̇yyp)
                # ∇vp  = -3/2*λ*dQdP

                # with this one: ezzp'!=0, ezzp=0, tzz!=0, ezz'=0
                ε̇zzp =  0.0
                ∇vp  = -λ*dQdP

                # ∇vp = ε̇xxp + ε̇yyp + ε̇zzp so ideally ε̇zzp = 0 such that ∇vp = -λ*dQdP

                ε̇xx = 0.
                ε̇yy = (0.333333333333333 * G .* Δt .* (-ε̇xxp + 2.0 * ε̇yyp - ε̇zzp + ∇vp) + 0.5 * P0 + 0.5 * σyyi - 0.5 * τyy0) ./ (G .* Δt)
                
                ε̇zz  = 1/2*(ε̇xx + ε̇yy)
                ∇v   = ε̇xx + ε̇yy + ε̇zz
                ε̇xxd = ε̇xx - 1/3*∇v
                ε̇yyd = ε̇yy - 1/3*∇v
                ε̇zzd = ε̇zz - 1/3*∇v

                τxx = τxx0 + 4/3*G*Δt*(ε̇xxd - ε̇xxp + 1/3*∇vp) - 2/3*G*Δt*(ε̇yyd - ε̇yyp + 1/3*∇vp) - 2/3*G*Δt*(ε̇zzd - ε̇zzp + 1/3*∇vp)
                τyy = τyy0 - 2/3*G*Δt*(ε̇xxd - ε̇xxp + 1/3*∇vp) + 4/3*G*Δt*(ε̇yyd - ε̇yyp + 1/3*∇vp) - 2/3*G*Δt*(ε̇zzd - ε̇zzp + 1/3*∇vp)
                τzz = τzz0 - 2/3*G*Δt*(ε̇xxd - ε̇xxp + 1/3*∇vp) - 2/3*G*Δt*(ε̇yyd - ε̇yyp + 1/3*∇vp) + 4/3*G*Δt*(ε̇zzd - ε̇zzp + 1/3*∇vp)
                τxy = τxy0 +     G*Δt*(γ̇xy  - γ̇xyp)
                
                P   = P0   - 2/3*G*Δt*(∇v- ∇vp) 
                
                ε̇xxd = ε̇xx - 1/3*(ε̇xx + ε̇yy + ε̇zz)
                ε̇yyd = ε̇yy - 1/3*(ε̇xx + ε̇yy + ε̇zz)
                ε̇zzd = ε̇zz - 1/3*(ε̇xx + ε̇yy + ε̇zz)

                σxx  = σxx0 + 2*G*Δt*(ε̇xx - ε̇xxp)
                σyy  = σyy0 + 2*G*Δt*(ε̇yy - ε̇yyp) # this will not change
                σxy  = σxy0 +   G*Δt*(γ̇xy - γ̇xyp)
                σzz  = σzz0 + 2*G*Δt*(ε̇zz - ε̇zzp)

                # σxx1 = σxx0 + 2*G*Δt*(ε̇xxd + 1/3*∇v - ε̇xxp)
                # σyy1 = σyy0 + 2*G*Δt*(ε̇yyd + 1/3*∇v - ε̇yyp)
                # σzz1 = σzz0 + 2*G*Δt*(ε̇zzd + 1/3*∇v - ε̇zzp)
                # P1   = -1/3*(σxx1+σyy1+σzz1)
                # P2   = P0 - 2/3*G*Δt*(ε̇xxd + 1/3*∇v - ε̇xxp) -2/3*G*Δt*(ε̇yyd + 1/3*∇v - ε̇yyp) -2/3*G*Δt*(ε̇zzd + 1/3*∇v - ε̇zzp)
                # P2   = P0 - 2/3*G*Δt*( ∇v - (ε̇xxp + ε̇yyp + ε̇zzp)  )
          
                τ    = sqrt(0.25*(τxx-τyy)^2 + τxy^2)
                fc   = τ + 0.5*(τxx+τyy-2P)*sin(ϕ) 
                λ   += fc / G / 10
                if abs(f)<1e-8 break end
            end
        end
        σyy = -P + τyy
        if it==30 @printf("%d %f %2.2e %f %2.2e %f %f %2.2e %2.2e %2.2e\n", it, f, λ, P, ε̇yy, σyy, τzz,   ε̇zz, ε̇zz-1/3*∇v, ∇v)             
            # @show  (-2P + τxx + τyy)/2, -P
            # @show -(-2P + τxx + τyy)/2, P
            # @show -1/3*(τxx+τyy+τzz-3P), P
            @printf("%1.2f %1.2f %1.2f %1.2f %1.2f\n", σxx, σyy, σzz, σxy, -1/3*(σxx+σyy+σzz))
            @show ε̇xx
            @show ε̇yy
            @show ε̇zz
            @show ε̇zzp - 1/3*∇vp
            @show τxx
            @show τyy
            @show τzz
            @show τzz-P
            @show τxx+τyy+τzz
            @show  ε̇zzd + ε̇xxd + ε̇yyd
        end
        # Storage
        app_fric[it] = -τxy/(τyy-P)
        γxy[it]      = γxy[it-1] + γ̇xy*Δt
        θv[it]       = θ
    end
    return γxy, θv, app_fric
end


function main()

    app_fric_V90_A = [0.3599745143682773 0.2952279659656446;
                    0.6731919248675551 0.42227484347407296;
                    1.0107812249157166 0.510816342912185;
                    1.5244170412586289 0.5892318189115429;
                    2.0122060122009957 0.6358604912506021;
                    2.9866250602022797 0.6738882645689517;
                    3.4988561566864664 0.6853588055867716;
                    4.011017017177718 0.6934820998555147;
                    4.4979631562048485 0.6999438112056511;
                    5.0015903435543425 0.7013766254615508;
                    5.4968243297479535 0.7028134532027613;
                    6.01302375983304 0.703403435543426;
                    6.4956855032910585 0.7056830951998716;
                    6.9992775726440835 0.705442286081233]


    app_fric_V90_B = [0.14708921977845546 0.14972306951356562;
                    0.3012923422700273 0.298603307111896;
                    0.5922800610049768 0.566246588537486;
                    1.0233885856477765 0.7116471343714882;
                    1.5204135896612614 0.7984387542141597;
                    2.007921616631883 0.8316784395569112;
                    2.5032960748113666 0.8398097607962756;
                    2.9982842350297005 0.829531224915717;
                    3.987979611494622 0.7955851661582919;
                    4.491220500882966 0.7786081232942688;
                    5.0029599454166 0.7666479370685504;
                    5.506235952801413 0.7513445175790657;
                    6.009652432172098 0.7427355915877348;
                    6.496282509231016 0.7341346925670253;
                    6.9997341065981695 0.7271993899502329;
                    7.486399301653557 0.7202721143040618;
                    7.9899562530101145 0.7183576818108846]

    params = (
        G   = 10e6,
        ϕ   = 40/180*π,
        ψ   = 10/180*π,
        γ̇xy = 0.00001,
        Δt  = 20,
        # nt  = 30,
        nt  = 400,
        pl  = true
    )

    style = :matrix
    # style = :componentwise

    #----------------------
    σxxA       = -25e3 # Courbe A - Vermeer
    σyyA       = -100e3 # Courbe A - Vermeer
    #----------------------
    σxxB       = -400e3 # Courbe B - Vermeer
    σyyB       = -100e3 # Courbe B - Vermeer
    #----------------------

    if style == :componentwise
        @info "Run with componentwise approach"
        γxy_A, θv_A, app_fric_A = Vermeer1990_FirstTest_Components(σxxA, σyyA, params)
        γxy_B, θv_B, app_fric_B = Vermeer1990_FirstTest_Components(σxxB, σyyB, params)
    else
        @info "Run with matrix-vector approach"
        γxy_A, θv_A, app_fric_A = Vermeer1990_FirstTest_Matrix(σxxA, σyyA, params)
        γxy_B, θv_B, app_fric_B = Vermeer1990_FirstTest_Matrix(σxxB, σyyB, params)
    end

    # Test your own stress integration scheme - total stress - Mohr-Coulonmb
    γxy_1, θv_1, app_fric_1 = Vermeer1990_StressIntegration1(σxxB, σyyB, params)

    # Test your own stress integration scheme - deviatoric stress / pressure - Mohr-Coulonmb
    γxy_2, θv_2, app_fric_2 = Vermeer1990_StressIntegration2(σxxB, σyyB, params)

    p1 = plot( xlabel = "γxy",  ylabel = "-σxy/σyy" )
    # p1 = plot!(γxy_A[2:end], app_fric_A[2:end], label="Case A - semi-analytics")
    p1 = plot!(γxy_B[2:end], app_fric_B[2:end], label="Case B - semi-analytics")
    p1 = plot!(γxy_1[5:6:end], app_fric_1[5:6:end], label="Numerics 1 MC", marker=:cross)
    p1 = scatter!(γxy_2[2:10:end], app_fric_2[2:10:end], label="Numerics 2 MC", marker=:star)
    
    p1 = scatter!(app_fric_V90_A[:,1]./100, app_fric_V90_A[:,2], label="Case A - Vermeer (1990)")
    p1 = scatter!(app_fric_V90_B[:,1]./100, app_fric_V90_B[:,2], label="Case B - Vermeer (1990)")



    p2 = plot( xlabel = "γxy",  ylabel = "θ (σ₃ angle)" )
    # p2 = plot!(γxy_A[2:end], θv_A[2:end]*180/π, label="case A - semi-analytics")
    p2 = plot!(γxy_B[2:end], θv_B[2:end]*180/π, label="case B - semi-analytics")
    p2 = scatter!(γxy_1[5:6:end], θv_1[5:6:end]*180/π, label="numerics 1 MC ", marker=:cross)
    p2 = scatter!(γxy_2[2:10:end], θv_2[2:10:end]*180/π, label="numerics 2 MC", marker=:star)

    display(plot(p1, p2, layout=(2,1)))
end

main()