# Define a mutable struct for parameters
mutable struct params
    # Simulation flags
    Etol::Float64
    rtol::Float64
    phi::Float64
    cut_off::Float64

    # Physical constants
    hbar::Float64
    kbol::Float64
    mu0::Float64
    muB::Float64
    a0::Float64
    m0::Float64
    w0::Float64
    mu0factor::Float64

    # System parameters
    m::Float64
    mu::Float64
    l0::Float64

    # Project info
    # name::String
    save_plots::Bool
    save_params::Bool
    save_intermediate_plots::Bool
    save_intermediate_data::Bool
    stop_relres_flag::Bool
    stop_relres::Float64
    # description::String

    # Grid / simulation
    Nx::Int
    Ny::Int
    Nz::Int
    Lx::Float64
    Ly::Float64
    Lz::Float64
    wx::Float64
    wy::Float64
    wz::Float64
    dt::Float64
    mindt::Float64

    # Variable parameters
    as_a0::Float64
    as::Float64
    T::Float64
    theta_deg::Float64
    theta::Float64
    N::Float64

    # Derived quantities
    gs::Float64
    add::Float64
    eps_dd::Float64
    Q5::Float64
    gammaQF::Float64
    gdd::Float64
    ellx::Float64
    elly::Float64
    ellz::Float64
    gx::Float64
    gy::Float64
    gz::Float64
end

# Factory function to create a Params object
function parameters()
    # Flags
    Etol = 5e-10
    rtol = 1e-6
    phi = 0.0
    cut_off = 3e7

    # Constants
    hbar = 1.0545718e-34
    kbol = 1.38064852e-23
    mu0  = 1.25663706212e-6
    muB  = 9.274009994e-24
    a0   = 5.2917721067e-11
    m0   = 1.660539066e-27
    w0   = 2*pi*61.63158647
    mu0factor = 0.3049584233607396

    # System
    m = 164*m0
    mu = 9.93*muB
    l0 = sqrt(hbar/(m*w0))

    # Project
    save_plots = true
    save_params = false
    save_intermediate_plots = true
    save_intermediate_data = true
    stop_relres_flag = false
    stop_relres = 10^(-4.5)
    # description = "Groundstate."

    # Grid
    Nx = 16
    Ny = 16
    Nz = 16
    Lx = 48
    Ly = 48
    Lz = 24
    wx = 2*pi*39
    wy = 2*pi*20
    wz = 2*pi*97
    dt = 0.001
    mindt = 2e-8

    # Variable parameters
    as_a0 = 97
    as = as_a0*a0
    T = 60
    theta_deg = 10
    theta = deg2rad(theta_deg)
    N = 200000

    # Derived quantities
    gs = 4*pi*as/l0
    add = mu0*mu^2*m/(12*pi*hbar^2)
    eps_dd = add/as
    Q5 = 1 + 3/2*eps_dd^2
    gammaQF = 128/3*sqrt(pi*(as/l0)^5)*Q5
    gdd = 4*pi*add/l0
    ellx = sqrt(hbar/(m*wx))
    elly = sqrt(hbar/(m*wy))
    ellz = sqrt(hbar/(m*wz))
    gx = (wx/w0)^2
    gy = (wy/w0)^2
    gz = (wz/w0)^2

    return params(
        Etol, rtol, phi, cut_off,
        hbar, kbol, mu0, muB, a0, m0, w0, mu0factor,
        m, mu, l0,
        save_plots, save_params, save_intermediate_plots, save_intermediate_data, stop_relres_flag, stop_relres, 
        Nx, Ny, Nz, Lx, Ly, Lz, wx, wy, wz, dt, mindt,
        as_a0, as, T, theta_deg, theta, N,
        gs, add, eps_dd, Q5, gammaQF, gdd, ellx, elly, ellz, gx, gy, gz
    )
end