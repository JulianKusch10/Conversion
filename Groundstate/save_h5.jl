using HDF5

"""
    save_h5(filename, name, x)

Saves `x` (array, number, or mutable/nested struct) to an HDF5 file `filename` 
under the group/dataset `name`. Works recursively for nested structs.
"""
function save_h5(filename::AbstractString, name::AbstractString, x)
    # inner recursive function
    function _write(parent, name, obj)
        if isstructtype(typeof(obj)) && !(obj isa AbstractArray)
            # mutable struct → save as a group
            g = create_group(parent, name)
            attrs(g)["julia_type"] = string(typeof(obj))
            for field in fieldnames(typeof(obj))
                _write(g, string(field), getfield(obj, field))
            end

        elseif eltype(obj) <: Complex
            # complex array → save real and imag separately
            parent["$(name)_real"] = real(obj)
            parent["$(name)_imag"] = imag(obj)

        else
            # normal numbers, arrays, strings
            parent[name] = obj
        end
    end

    # open the file: "r+" if exists, "w" if not
    mode = isfile(filename) ? "r+" : "w"
    h5open(filename, mode) do f
        _write(f, name, x)
    end
end