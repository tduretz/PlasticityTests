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