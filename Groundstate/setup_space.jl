# setup_space.jl
using FFTW

mutable struct TransfType
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    kx::Vector{Float64}
    ky::Vector{Float64}
    kz::Vector{Float64}
    X::Array{Float64,3}
    Y::Array{Float64,3}
    Z::Array{Float64,3}
    KX::Array{Float64,3}
    KY::Array{Float64,3}
    KZ::Array{Float64,3}
    dx::Float64
    dy::Float64
    dz::Float64
    dkx::Float64
    dky::Float64
    dkz::Float64
end

function setup_space(Params)
    Nx, Ny, Nz = Params.Nx, Params.Ny, Params.Nz

    # Real space grids
    x = range(-0.5*Params.Lx, 0.5*Params.Lx, length=Nx)
    y = range(-0.5*Params.Ly, 0.5*Params.Ly, length=Ny)
    z = range(-0.5*Params.Lz, 0.5*Params.Lz, length=Nz)

    dx = x[2]-x[1]
    dy = y[2]-y[1]
    dz = z[2]-z[1]

    # Fourier grids
    kx = fftshift(2*pi*fftfreq(Nx, dx))
    ky = fftshift(2*pi*fftfreq(Ny, dy))
    kz = fftshift(2*pi*fftfreq(Nz, dz))

    dkx = kx[2]-kx[1]
    dky = ky[2]-ky[1]
    dkz = kz[2]-kz[1]

    # Meshgrids
    X = reshape(x, Nx, 1, 1) .* ones(1, Ny, Nz)
    Y = reshape(y, 1, Ny, 1) .* ones(Nx, 1, Nz)
    Z = reshape(z, 1, 1, Nz) .* ones(Nx, Ny, 1)

    KX = reshape(kx, Nx, 1, 1) .* ones(1, Ny, Nz)
    KY = reshape(ky, 1, Ny, 1) .* ones(Nx, 1, Nz)
    KZ = reshape(kz, 1, 1, Nz) .* ones(Nx, Ny, 1)

    return TransfType(x, y, z, kx, ky, kz, X, Y, Z, KX, KY, KZ, dx, dy, dz, dkx, dky, dkz)
end