function print_structures(p::ParamsType, filename::AbstractString)
    open(filename, "w") do io
        for field in fieldnames(typeof(p))
            value = getfield(p, field)
            if isa(value, AbstractArray)
                println(io, string(field, " = array of size ", size(value)))
            else
                println(io, string(field, " = ", value))
            end
        end
    end
end