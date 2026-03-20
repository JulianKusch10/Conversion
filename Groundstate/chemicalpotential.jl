function chemicalpotential(psi, Params, Transf, VDk, V)

    # Parameters
    normfac = Params.Lx * Params.Ly * Params.Lz / length(psi)
    KEop = 0.5 .* (Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)

    # DDI
    frho = fft(abs.(psi).^2)
    Phi = ifft(frho .* VDk)

    Eddi = Params.gdd .* Phi .* abs.(psi).^2
    Eddi = sum(Eddi) * Transf.dx * Transf.dy * Transf.dz

    # Kinetic energy
    Ekin = conj.(psi) .* ifft(KEop .* fft(psi))
    Ekin = sum(Ekin) * Transf.dx * Transf.dy * Transf.dz

    # Potential energy
    Epot = V .* abs.(psi).^2
    Epot = sum(Epot) * Transf.dx * Transf.dy * Transf.dz

    # Contact interactions
    Eint = Params.gs .* abs.(psi).^4
    Eint = sum(Eint) * Transf.dx * Transf.dy * Transf.dz

    # Quantum fluctuations
    Eqf = Params.gammaQF .* abs.(psi).^5
    Eqf = sum(Eqf) * Transf.dx * Transf.dy * Transf.dz

    # Thermal energy
    Eth = 0
    # Eth = HT .* abs.(psi).^2
    # Eth = sum(Eth) * Transf.dx * Transf.dy * Transf.dz

    # Chemical potential
    muchem = real(Ekin + Epot + Eint + Eddi + Eqf + Eth) / Params.N

    return muchem
end