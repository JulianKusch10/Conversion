using HDF5
include("parameters.jl")
include("setup_space.jl")

"""
    load_h5(filename, name)

Loads data from an HDF5 file saved with `save_h5`.  
Automatically reconstructs mutable structs if the type exists.
"""
function load_h5(filename::AbstractString, name::AbstractString)

    function _read(obj)
        if obj isa HDF5.Group
            # check if this group has a stored Julia type
            T = haskey(attrs(obj), "julia_type") ? eval(Meta.parse(attrs(obj)["julia_type"])) : nothing

            # recursively read all fields
            vals = Dict{String,Any}()
            for k in keys(obj)
                vals[string(k)] = _read(obj[k])
            end

            if T !== nothing
                # construct struct using field order
                args = [vals[string(f)] for f in fieldnames(T)]
                return T(args...)
            else
                return vals  # fallback to Dict if no type info
            end
        else
            return read(obj)  # dataset
        end
    end

    h5open(filename, "r") do f
        return _read(f[name])
    end
end