using Plots

function BDF_VE_relaxation(bdf, n)

    η  = 10.0
    G  = 1.0
    Δt = 1.0/n
    Nt = 100*n
    ε̇  = 1.0

    t0 = 0.
    τ     = 2η*ε̇*(1 .- exp.(-G/η.*t0))
    τ0    = 2η*ε̇*(1 .- exp.(-G/η.*(t0-1*Δt)))
    τ00   = 2η*ε̇*(1 .- exp.(-G/η.*(t0-2*Δt)))
    τ000  = 2η*ε̇*(1 .- exp.(-G/η.*(t0-3*Δt)))
    τ0000 = 2η*ε̇*(1 .- exp.(-G/η.*(t0-4*Δt)))

    # τ     = 2η*ε̇*(1 .- exp.(-G/η.*t0))
    # τ0    = 2η*ε̇*(1 .- exp.(-G/η.*(t0-0*1*Δt)))
    # τ00   = 2η*ε̇*(1 .- exp.(-G/η.*(t0-0*2*Δt)))
    # τ000  = 2η*ε̇*(1 .- exp.(-G/η.*(t0-0*3*Δt)))
    # τ0000 = 2η*ε̇*(1 .- exp.(-G/η.*(t0-0*4*Δt)))

    τvec = zeros(Nt)

    a, b, c, d, e = 1/Δt, -1/Δt, 0., 0., 0.
    
    for it=1:Nt

        if bdf==2
            a = 3/2/Δt
            b =  -2/Δt
            c = 1/2/Δt
            d =   0/Δt
            e =   0/Δt
        end
        if bdf==3
            a = 11/6/Δt
            b =   -3/Δt
            c =  3/2/Δt
            d = -1/3/Δt
            e =   0/Δt
        end
        if bdf==4
            a =  25/12/Δt
            b = -48/25*a
            c =  36/25*a
            d = -16/25*a
            e =   3/25*a
        end

        # if it>1 && bdf>1
        #     a = 3/2/Δt
        #     b =  -2/Δt
        #     c = 1/2/Δt
        #     d =   0/Δt
        #     e =   0/Δt
        # end
        # if it>2 && bdf>2
        #     a = 11/6/Δt
        #     b =   -3/Δt
        #     c =  3/2/Δt
        #     d = -1/3/Δt
        #     e =   0/Δt
        # end
        # if it>3 && bdf>3
        #     a =  25/12/Δt
        #     b = -48/25*a
        #     c =  36/25*a
        #     d = -16/25*a
        #     e =   3/25*a
        # end

        ηve  = ( 1/η + a/G )^-1

        τ0000 = τ000
        τ000  = τ00
        τ00   = τ0
        τ0    = τ
        τ     = 2*ηve* (ε̇ - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        τvec[it] = τ
    end

    p   = plot()
    t   = (1:Nt).*Δt
    τvec_ana = 2η*ε̇*(1 .- exp.(-G/η.*t))
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

@info "BDF2"
res[2,1] = BDF_VE_relaxation(2,1)
res[2,2] = BDF_VE_relaxation(2,2)
res[2,3] = BDF_VE_relaxation(2,4)
res[2,4] = BDF_VE_relaxation(2,8)

@info "BDF3"
res[3,1] = BDF_VE_relaxation(3,1)
res[3,2] = BDF_VE_relaxation(3,2)
res[3,3] = BDF_VE_relaxation(3,4)
res[3,4] = BDF_VE_relaxation(3,8)

@info "BDF4"
res[4,1] = BDF_VE_relaxation(4,1)
res[4,2] = BDF_VE_relaxation(4,2)
res[4,3] = BDF_VE_relaxation(4,4)
res[4,4] = BDF_VE_relaxation(4,8)
@show res[4,:]

x = 1 ./[1 2 4 8]
p=plot()
p=plot!(log10.(1 ./x[:]), log10.(res[1,:]), marker=:dot, label="BDF1")
p=plot!(log10.(1 ./x[:]), log10.(res[2,:]), marker=:dot, label="BDF2")
p=plot!(log10.(1 ./x[:]), log10.(res[3,:]), marker=:dot, label="BDF3")
p=plot!(log10.(1 ./x[:]), log10.(res[4,:]), marker=:dot, label="BDF4")
