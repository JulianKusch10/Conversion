Base.@kwdef mutable struct ObservablesType
    # Example observables
    EVec::Vector{Float64} = Float64[]
    mucVec::Vector{Float64} = Float64[]
    residual::Vector{Float64} = Float64[]
    # You can add more fields as needed
end

# Factory function to create an initialized observables object
function observs()
    return ObservablesType()  # returns object with default values
end