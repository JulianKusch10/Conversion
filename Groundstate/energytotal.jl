function energytotal(psi, Params, Transf, VDk, V, Ftherm)

    # Parameters
    KEop = 0.5 .* (Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)

    # Dipole-dipole interactions
    frho = fft(abs.(psi).^2)
    Phi = ifft(frho .* VDk)

    Eddi = 0.5 * Params.gdd .* Phi .* abs.(psi).^2
    Eddi = sum(Eddi) * Transf.dx * Transf.dy * Transf.dz

    # Kinetic energy
    Ekin = conj.(psi) .* ifft(KEop .* fft(psi))
    Ekin = sum(Ekin) * Transf.dx * Transf.dy * Transf.dz

    # Potential energy
    Epot = V .* abs.(psi).^2
    Epot = sum(Epot) * Transf.dx * Transf.dy * Transf.dz

    # Contact interactions
    Eint = 0.5 * Params.gs .* abs.(psi).^4
    Eint = sum(Eint) * Transf.dx * Transf.dy * Transf.dz

    # Quantum fluctuations
    Eqf = 0.4 * Params.gammaQF .* abs.(psi).^5
    Eqf = sum(Eqf) * Transf.dx * Transf.dy * Transf.dz

    # Thermal energy
    n = abs.(psi).^2
    Fth = Ftherm.(n)
    Eth = sum(Fth) * Transf.dx * Transf.dy * Transf.dz

    return real(Ekin + Epot + Eint + Eddi + Eqf + Eth)

end