@doc raw"""
    Case = ExtractDataCase(CaseName)  ;

Function that reads all CSV files containing the digitized data from the test 1 in Vermeer (1990). 
The function takes as inputs:

    CaseName : Name of the case, string ("CaseA" or "CaseB") 
and returns:

    Case : A tuple that contains all digitized data

# Examples
```julia-repl
julia>  CaseA = ExtractDataCase("CaseA") 

```
"""
@views function ExtractDataCase( CaseName )
    Case = (
        Friction = (
            x = convert(Array, CSV.read("./data/Vermeer1990_Friction_$(CaseName).csv", DataFrame, header =[:x, :y])[:,1]),
            y = convert(Array, CSV.read("./data/Vermeer1990_Friction_$(CaseName).csv", DataFrame, header =[:x, :y])[:,2]),
        ),
        σxx = (
            x = convert(Array, CSV.read("./data/Vermeer1990_Sxx_$(CaseName).csv", DataFrame, header =[:x, :y])[:,1]),
            y = convert(Array, CSV.read("./data/Vermeer1990_Sxx_$(CaseName).csv", DataFrame, header =[:x, :y])[:,2]),
        ),
        εyy = (
            x = convert(Array, CSV.read("./data/Vermeer1990_Eyy_$(CaseName).csv", DataFrame, header =[:x, :y])[:,1]),
            y = convert(Array, CSV.read("./data/Vermeer1990_Eyy_$(CaseName).csv", DataFrame, header =[:x, :y])[:,2]),
        ),
        θ = (
            x = convert(Array, CSV.read("./data/Vermeer1990_theta_$(CaseName).csv", DataFrame, header =[:x, :y])[:,1]),
            y = convert(Array, CSV.read("./data/Vermeer1990_theta_$(CaseName).csv", DataFrame, header =[:x, :y])[:,2]),
        ),
    )
    return Case
end