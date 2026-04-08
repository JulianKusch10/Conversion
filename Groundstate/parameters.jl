############### parameters.jl ################
Base.@kwdef mutable struct ParamsType

    # --------------------------
    # Simulation flags
    # --------------------------
    Etol::Float64 = 5e-10
    rtol::Float64 = 1e-6
    phi::Float64 = 0.0
    cut_off::Float64 = 3e7

    # --------------------------
    # Physical constants
    # --------------------------
    hbar::Float64 = 1.0545718e-34
    kbol::Float64 = 1.38064852e-23
    mu0::Float64 = 1.25663706212e-6
    muB::Float64 = 9.274009994e-24
    a0::Float64 = 5.2917721067e-11
    m0::Float64 = 1.660539066e-27
    w0::Float64 = 2π*61.63158647
    mu0factor::Float64 = 0.3049584233607396

    # --------------------------
    # System parameters
    # --------------------------
    m::Float64 = 164*m0
    mu::Float64 = 9.93*muB
    l0::Float64 = sqrt(hbar/(m*w0))

    # --------------------------
    # Project settings
    # --------------------------
    save_plots::Bool = true
    save_params::Bool = false
    save_intermediate_plots::Bool = true
    save_intermediate_data::Bool = true
    stop_relres_flag::Bool = false
    stop_relres::Float64 = 10^(-4.5)

    # --------------------------
    # Grid / simulation
    # --------------------------
    Nx::Int = 16
    Ny::Int = 16
    Nz::Int = 16

    Lx::Float64 = 48
    Ly::Float64 = 48
    Lz::Float64 = 24

    wx::Float64 = 2π*39
    wy::Float64 = 2π*20
    wz::Float64 = 2π*97

    dt::Float64 = 0.001
    mindt::Float64 = 2e-8

    # --------------------------
    # Variable parameters
    # --------------------------
    as_a0::Float64 = 97
    as::Float64 = as_a0*a0
    T::Float64 = 60
    theta_deg::Float64 = 10
    theta::Float64 = deg2rad(theta_deg)
    N::Float64 = 200000

    # --------------------------
    # Derived quantities
    # --------------------------
    gs::Float64 = NaN
    add::Float64 = NaN
    eps_dd::Float64 = NaN
    Q5::Float64 = NaN
    gammaQF::Float64 = NaN
    gdd::Float64 = NaN

    ellx::Float64 = NaN
    elly::Float64 = NaN
    ellz::Float64 = NaN

    gx::Float64 = NaN
    gy::Float64 = NaN
    gz::Float64 = NaN

    fx::Float64 = NaN
    fy::Float64 = NaN
    fz::Float64 = NaN

end


# ============================================================
# Function to update derived quantities
# ============================================================

function update_params!(Params::ParamsType)

    # variable conversions
    Params.as = Params.as_a0 * Params.a0
    Params.theta = deg2rad(Params.theta_deg)

    # derived physics
    Params.gs = 4π*Params.as/Params.l0
    Params.add = Params.mu0*Params.mu^2*Params.m/(12π*Params.hbar^2)
    Params.eps_dd = Params.add/Params.as
    Params.Q5 = 1 + 3/2*Params.eps_dd^2
    Params.gammaQF = 128/3*sqrt(π*(Params.as/Params.l0)^5)*Params.Q5
    Params.gdd = 4π*Params.add/Params.l0

    # oscillator lengths
    Params.ellx = sqrt(Params.hbar/(Params.m*Params.wx))
    Params.elly = sqrt(Params.hbar/(Params.m*Params.wy))
    Params.ellz = sqrt(Params.hbar/(Params.m*Params.wz))

    # trap strengths
    Params.gx = (Params.wx/Params.w0)^2
    Params.gy = (Params.wy/Params.w0)^2
    Params.gz = (Params.wz/Params.w0)^2

    # Conversion to frequency
    Params.fx = Params.wx/(2*pi)
    Params.fy = Params.wy/(2*pi)
    Params.fz = Params.wz/(2*pi)

end


# ============================================================
# Factory function returning Params object
# ============================================================

function parameters()
    Params = ParamsType()
    update_params!(Params)
    return Params
end