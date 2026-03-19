using LinearAlgebra
using FFTW
using Random
using Statistics
using Interpolations
using SpecialFunctions

function Initialize(Params, Transf)

    # == Grid ==
    X, Y, Z = Transf.X, Transf.Y, Transf.Z

    # == Potential ==
    V = 0.5 .* (Params.gx .* X.^2 .+ Params.gy .* Y.^2 .+ Params.gz .* Z.^2)

    # == Calculating the DDIs ==
    CutoffType = 6
    VDk = nothing

    if CutoffType == 1
        # Cylindrical (semi-analytic)
        Zcutoff = Params.Lz / 2
        alph = acos.(
            (Transf.KX .* sin(Params.theta) * cos(Params.phi) .+ 
             Transf.KY .* sin(Params.theta) * sin(Params.phi) .+ 
             Transf.KZ .* cos(Params.theta)) ./ 
            sqrt.(Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)
        )
        alph[1] = π/2

        # Analytic part
        cossq = cos.(alph).^2
        VDk = cossq .- 1/3
        sinsq = 1 .- cossq
        VDk .+= exp.(-Zcutoff .* sqrt.(Transf.KX.^2 .+ Transf.KY.^2)) .* 
               (sinsq .* cos.(Zcutoff .* Transf.KZ) .- sqrt.(sinsq .* cossq) .* sin.(Zcutoff .* Transf.KZ))

        # Nonanalytic part
        Params.Lr = 0.5 * min(Params.Lx, Params.Ly)
        Params.Nr = max(Params.Nx, Params.Ny)
        TransfRad = setup_space_radial(Params)
        VDkNon = VD_RadInt(TransfRad.kr, TransfRad.kz, TransfRad.Rmax, Zcutoff)

        # Interpolation onto 3D grid
        fullkr = vcat(reverse(TransfRad.kr), TransfRad.kr)
        KR, KZ = TransfRad.kr, TransfRad.kz
        # 3D grid (KX3D,KY3D,KZ3D) via broadcasting
        KR3D = sqrt.(ifftshift.(Transf.kx).^2 .+ ifftshift.(Transf.ky).^2)
        fullVDK = hcat(reverse(VDkNon, dims=2), VDkNon)
        # Interpolation (using Interpolations.jl)
        itp = interpolate((KR, KZ), fullVDK, Gridded(Linear()))
        VDkNon = itp.(KR3D, KZ)
        VDkNon = fftshift(VDkNon)

        VDk .+= VDkNon
        VDk .*= 3  # Fixing 1/3 inconsistency

        # Optionally save
        # JLD2.jl or MAT.jl can be used to save as .mat
        println("Finished DDI")

    elseif CutoffType == 2
        # Cylindrical infinite Z, polarized along x
        alph = acos.(
            (Transf.KX .* sin(Params.theta) * cos(Params.phi) .+ 
             Transf.KY .* sin(Params.theta) * sin(Params.phi) .+ 
             Transf.KZ .* cos(Params.theta)) ./ 
            sqrt.(Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)
        )
        alph[1] = π/2
        rhoc = maximum(abs.(vcat(Transf.x, Transf.y)))
        KR = sqrt.(Transf.KX.^2 .+ Transf.KY.^2)

        func(n,u,v) = (v.^2 ./ (u.^2 .+ v.^2)) .* (v .* besselj.(n,u) .* besselk.(n+1,v) .- u .* besselj.(n+1,u) .* besselk.(n,v))
        VDk = -0.5 .* func(0, KR .* rhoc, abs.(Transf.KZ) .* rhoc) .+ 
              (Transf.KX.^2 ./ KR.^2 .- 0.5) .* func(2, KR .* rhoc, abs.(Transf.KZ) .* rhoc)
        VDk = (1/3) .* (3 .* cos.(alph).^2 .- 1) .- VDk

        VDk[KR .== 0] .= -1/3 .+ 0.5 .* abs.(Transf.KZ[KR .== 0]) .* rhoc .* besselk.(1, abs.(Transf.KZ[KR .== 0]) .* rhoc)
        VDk[Transf.KZ .== 0] .= 1/6 .+ (Transf.KX[Transf.KZ .== 0].^2 .- Transf.KY[Transf.KZ .== 0].^2) ./ (KR[Transf.KZ .== 0].^2) .* (1/2 .- besselj.(1, KR[Transf.KZ .== 0] .* rhoc) ./ (KR[Transf.KZ .== 0] .* rhoc))
        VDk[1,1,1] = 1/6
        VDk .*= 3

    elseif CutoffType == 6
        # Rectangular
        alph = acos.(
            (Transf.KX .* sin(Params.theta) * cos(Params.phi) .+ 
             Transf.KY .* sin(Params.theta) * sin(Params.phi) .+ 
             Transf.KZ .* cos(Params.theta)) ./ 
            sqrt.(Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)
        )
        alph[1] = π/2
        VDk = cos.(alph).^2 .- 1/3
        VDk[1,1,1] = 0

        Xcutoff, Ycutoff, Zcutoff = Params.Lx/2, Params.Ly/2, Params.Lz/2
        VDr = fftshift(fft(VDk))
        VDr[abs.(Transf.X) .> Xcutoff] .= 0
        VDr[abs.(Transf.Y) .> Ycutoff] .= 0
        VDr[abs.(Transf.Z) .> Zcutoff] .= 0
        VDk = ifft(ifftshift(VDr))
        VDk .*= 3

    else
        println("Choose a valid DDI!")
        return nothing
    end

    # == Initial wavefunction ==
    
        ellx = sqrt(Params.hbar / (Params.m * Params.wx)) / Params.l0
        elly = sqrt(Params.hbar / (Params.m * Params.wy)) / Params.l0
        ellz = sqrt(Params.hbar / (Params.m * Params.wz)) / Params.l0

        Rx, Ry, Rz = 8*ellx, 8*elly, 8*ellz
        X0, Y0, Z0 = 0, 0, 0

        psi = exp.(-(X.-X0).^2 ./ Rx^2 .- (Y.-Y0).^2 ./ Ry^2 .- (Z.-Z0).^2 ./ Rz^2)
        cur_norm = sum(abs2.(psi)) * Transf.dx * Transf.dy * Transf.dz
        psi ./= sqrt(cur_norm)

        # Add noise
        # r = randn(size(X))
        # theta = rand(size(X))
        # noise = r
        # psi .+= 0.005 .* noise

        # Norm = sum(abs2.(psi)) * Transf.dx * Transf.dy * Transf.dz
        # psi .= sqrt(Params.N) .* psi ./ sqrt(Norm)

    return psi, V, VDk
end