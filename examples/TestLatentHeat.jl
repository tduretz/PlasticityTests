using Plots

ϕf(t,τk)    = 1 - exp(-t/τk)

function main()
    τk  = 1e9
    nt  = 100
    t   = 0 
    ϕ   = ϕf.(t, τk) 
    cp1 = 1000.
    cp2 = 1050.
    ρ1  = 3000.
    ρ2  = 3200.
    T   =  773.
    U   = ϕ*cp1*ρ1*T + (1-ϕ)*cp2*ρ2*T
    Δt  = 1e8

    hist = (t=zeros(nt), T=zeros(nt), ϕ=zeros(nt))
    
    for it=1:nt
        t  = t + Δt
        U0 = U
        for iter=1:50
            ϕ  = ϕf.(t, τk) 
            U  = ϕ*cp1*ρ1*T + (1-ϕ)*cp2*ρ2*T
            R  = (U - U0)/Δt
            T -= R*50
        end
        hist.t[it] = t
        hist.T[it] = T
        hist.ϕ[it] = ϕ
    end
    plot(hist.t, hist.T)
end

main()