using Plots

function bdf2(dt0, dt)
    a = 1 ./ (dt + dt0) + 1 ./ dt
    b = dt .* (-1 ./ dt0 - 1 ./ dt) ./ (dt + dt0) - 1 ./ dt
    c = dt ./ (dt0 .* (dt + dt0))
    # a= 3/2/dt
    # b = -2/dt
    # c = 1/2/dt
    return (a, b, c)
end

function bdf3(dt00, dt0, dt)
    a = 1 ./ (dt + dt0 + dt00) + 1 ./ (dt + dt0) + 1 ./ dt
    b = dt .* (dt + dt0) .* ((-1 ./ dt0 - 1 ./ dt) ./ (dt + dt0) - 1 ./ (dt0 .* (dt0 + dt00))) ./ (dt + dt0 + dt00) + dt .* (-1 ./ dt0 - 1 ./ dt) ./ (dt + dt0) - 1 ./ dt
    c = dt .* (dt + dt0) .* (-(-1 ./ dt00 - 1 ./ dt0) ./ (dt0 + dt00) + 1 ./ (dt0 .* (dt + dt0))) ./ (dt + dt0 + dt00) + dt ./ (dt0 .* (dt + dt0))
    d = -dt .* (dt + dt0) ./ (dt00 .* (dt0 + dt00) .* (dt + dt0 + dt00))
    return (a, b, c, d)
end

function BDF_VE_relaxation(bdf, n)
    
    η      = 10.0
    G      = 1.0
    Δt_ref = 1.0/n
    Nt     = 100*n
    ε̇      = 1.0
    C      = 10

    t0 = 0.
    τ     = 2η*ε̇*(1 .- exp.(-G/η.*t0))
    τ0    = 2η*ε̇*(1 .- exp.(-G/η.*(t0-1*Δt_ref)))
    τ00   = 2η*ε̇*(1 .- exp.(-G/η.*(t0-2*Δt_ref)))
    τ000  = 2η*ε̇*(1 .- exp.(-G/η.*(t0-3*Δt_ref)))
    τ0000 = 2η*ε̇*(1 .- exp.(-G/η.*(t0-4*Δt_ref)))

    τvec = zeros(Nt)

    # MAKE 1 TINY step
    Δt = Δt_ref/1e7
    a, b, c, d, e = 1/Δt, -1/Δt, 0., 0., 0.
    Δt0 = Δt

    for it=1:2
        τ00   = τ0
        τ0    = τ       
        ηve   = ( 1/η + a/G )^-1
        τ     = 2*ηve* (ε̇ - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        τvec[it] = τ
    end

  
    for it=1:Nt

        τ0000 = τ000
        τ000  = τ00
        τ00   = τ0
        τ0    = τ

        Δt00  = Δt0
        Δt0   = Δt
        Δt    = Δt_ref
        a, b, c, d, e = 1/Δt, -1/Δt, 0., 0., 0.

        if bdf==2
            a, b, c = bdf2(Δt0, Δt)
            d =   0/Δt
            e =   0/Δt
        end
        if bdf==3
            a, b, c, d = bdf3(Δt00, Δt0, Δt)
            # a = 11/6/Δt
            # b =   -3/Δt
            # c =  3/2/Δt
            # d = -1/3/Δt
            # e =   0/Δt
        end
        if bdf==4
            a =  25/12/Δt
            b = -48/25*a
            c =  36/25*a
            d = -16/25*a
            e =   3/25*a
        end
        ηve   = ( 1/η + a/G )^-1
        τ     = 2*ηve* (ε̇ - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        F     = abs(τ) - C
        if F>0
            λ̇     = 0.5*F/ηve
            ε̇pl   = λ̇
            τ     = 2*ηve* (ε̇ - ε̇pl - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        end
        τvec[it] = τ
    end

    t   = (1:Nt).*Δt
    τvec_ana = 2η*ε̇*(1 .- exp.(-G/η.*t))
    τvec_ana[τvec_ana.>C] .= 10
    
    p = plot()
    p = plot!(t, τvec)
    p = plot!(t, τvec_ana )
    display(p)

    err = norm(τvec .- τvec_ana)/norm(τvec_ana)
    return err

end

res = zeros(4,4)

@info "BDF1"
res[1,1] = BDF_VE_relaxation(1,1)
res[1,2] = BDF_VE_relaxation(1,2)
res[1,3] = BDF_VE_relaxation(1,4)
res[1,4] = BDF_VE_relaxation(1,8)
@show res[1,:]

@info "BDF2"
res[2,1] = BDF_VE_relaxation(2,1)
res[2,2] = BDF_VE_relaxation(2,2)
res[2,3] = BDF_VE_relaxation(2,4)
res[2,4] = BDF_VE_relaxation(2,8)
@show res[2,:]

@info "BDF3"
res[3,1] = BDF_VE_relaxation(3,1)
res[3,2] = BDF_VE_relaxation(3,2)
res[3,3] = BDF_VE_relaxation(3,4)
res[3,4] = BDF_VE_relaxation(3,8)
@show res[3,:]

# @info "BDF4"
# res[4,1] = BDF_VE_relaxation(4,1)
# res[4,2] = BDF_VE_relaxation(4,2)
# res[4,3] = BDF_VE_relaxation(4,4)
# res[4,4] = BDF_VE_relaxation(4,8)
# @show res[4,:]

x = 1 ./[1 2 4 8]
p=plot()
p=plot!(log10.(1 ./x[:]), log10.(res[1,:]), marker=:dot, label="BDF1")
p=plot!(log10.(1 ./x[:]), log10.(res[2,:]), marker=:dot, label="BDF2")
p=plot!(log10.(1 ./x[:]), log10.(res[3,:]), marker=:dot, label="BDF3")
# p=plot!(log10.(1 ./x[:]), log10.(res[4,:]), marker=:dot, label="BDF4")
