using MAT

include("load_h5.jl")
mat = matread("C:/data/groundstates_julia_T4_200000_60nK_01_60_200000_0_94.0.mat")
# Params.as_a0
# Params_loaded = load_h5("./compdata/groundstate.h5", "Params");
# Transf_loaded = load_h5("./compdata/groundstate.h5", "Transf");
# psi_real_loaded = load_h5("./compdata/groundstate.h5", "psi_real");
# psi_imag_loaded = load_h5("./compdata/groundstate.h5", "psi_imag");

# f = h5open("./compdata/groundstate.h5", "r")  # read-only