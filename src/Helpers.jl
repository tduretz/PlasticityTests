_minmax(name, var::Array{Int64})   = @sprintf("min(%06s) = %+011d --- max(%06s) = %+011d\n", name, minimum(var), name, maximum(var))
_minmax(name, var::Array{Float64}) = @sprintf("min(%06s) = %11.3e --- max(%06s) = %11.3e\n", name, minimum(var), name, maximum(var))

@doc raw"""
    @minmax(x)  

Macro that prints minimum and maximum value of an array:

    x   : array of Float64 or Int64 
and prints th minimum and maximum value of the array.

# Examples
```julia-repl
julia>  @minmax(x) 

```
"""
macro minmax(x)
    quote
        name  = $(string(x))
        var   = $(esc(x))
        info  = _minmax(name, var)
        print(info)
    end
end


function PrincipalStress!(σ1, σ3, τxx, τyy, τzz, τxy, P)
    for i in eachindex(τxy)
        σ  = @SMatrix[-P[i]+τxx[i] τxy[i] 0.; τxy[i] -P[i]+τyy[i] 0.; 0. 0. -P[i]+τzz[i]]
        v  = eigvecs(σ)
        σp = eigvals(σ)
        σ1.x[i] = v[1,1]
        σ1.z[i] = v[2,1]
        σ3.x[i] = v[1,3]
        σ3.z[i] = v[2,3]
        σ1.v[i] = σp[1]
        σ3.v[i] = σp[3]
    end
end