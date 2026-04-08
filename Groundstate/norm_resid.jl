using FFTW
function norm_resid(psi, Params, Transf, VDk, V, muchem, HT)

    KEop = 0.5 .* (Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)

    frho = fft(abs.(psi).^2)
    Phi  = ifft(frho .* VDk)

    Eddi = Params.gdd .* Phi .* psi
    Ekin = ifft(KEop .* fft(psi))
    Epot = V .* psi
    Eint = Params.gs .* abs.(psi).^2 .* psi
    Eqf  = Params.gammaQF .* abs.(psi).^3 .* psi
    Eth  = HT .* psi

    dV = Transf.dx * Transf.dy * Transf.dz

    num = sum(abs.(Ekin .+ Epot .+ Eint .+ Eddi .+ Eqf .+ Eth .- muchem .* psi)) * dV
    den = sum(abs.(muchem .* psi)) * dV

    return num / den
end