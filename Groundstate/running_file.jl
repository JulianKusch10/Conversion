include("parameters.jl")
include("setup_space.jl")
include("initialize.jl")
include("ssfm_imag.jl")
include("HTherm.jl")
include("EthermInterp.jl")
include("energytotal.jl")
include("norm_resid.jl")
include("runningplot.jl")
include("observables.jl")
include("chemicalpotential.jl")
include("save_h5.jl")
include("print_structures.jl")
using HDF5 
t_start = time()
println("---------")
println("Starting")


Params = parameters();          # create mutable parameters
Transf = setup_space(Params);   # create grids
psi, V, VDk = Initialize(Params, Transf);
Observ = observs()

h5open("./compdata/initialize_julia.h5", "w") do file
    write(file, "psi", psi)
    write(file, "V", V)
    write(file, "VDk", VDk)
    write(file, "dx", Transf.dx)
    write(file, "dy", Transf.dy)
    write(file, "dz", Transf.dz)
end

t_idx = 1

ssfm_imag(psi, Params, Transf, VDk, V, t_idx, Observ)


println("Finished")
println("---------")
t_elapsed = time() - t_start
println("Elapsed time: ", t_elapsed, " seconds")