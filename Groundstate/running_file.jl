include("parameters.jl")
include("setup_space.jl")
include("initialize.jl")
using HDF5 

Params = parameters();          # create mutable parameters
Transf = setup_space(Params);   # create grids
psi, V, VDk = Initialize(Params, Transf);

h5open("./compdata/initialize_julia.h5", "w") do file
    write(file, "psi", psi)
    write(file, "V", V)
    write(file, "VDk", VDk)
    write(file, "dx", Transf.dx)
    write(file, "dy", Transf.dy)
    write(file, "dz", Transf.dz)
end