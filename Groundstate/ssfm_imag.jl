function ssfm_imag(psi, Params, Transf, VDk, V, t_idx)

    dt = -1im * abs(Params.dt)

    KEop = 0.5 .* (Transf.KX.^2 .+ Transf.KY.^2 .+ Transf.KZ.^2)

    # Observ.residual = [1e-2]
    # Observ.res = 1
    # Observ.EddiVec = [10.0]

    HT = HTherm(Params, psi)
    muchem = chemicalpotential(psi, Params, Transf, VDk, V, HT)
    # runningplot(psi, Params, Transf, Observ)

    # Ftherm = EthermInterp(Params, HThermFit)

    # AdaptIdx = 0
    # psi = real(psi)

    while t_idx < 100 #Params.cut_off

        # Kinetic half-step
        psi = fft(psi)
        psi = psi .* exp.(-0.5im .* dt .* KEop)
        psi = ifft(psi)

        # DDI
        frho = fft(abs.(psi).^2)
        Phi = ifft(frho .* VDk)

        # Real-space step
        # HT = HTherm(HThermFit, Params, psi)
        HT = HTherm(Params, psi)
        psi = psi .* exp.(-1im .* dt .* (V .+ Params.gs .* abs.(psi).^2 .+
                                          Params.gammaQF .* abs.(psi).^3 .+
                                          Params.gdd .* Phi .+ HT .- muchem))

        # Kinetic half-step
        psi = fft(psi)
        psi = psi .* exp.(-0.5im .* dt .* KEop)
        psi = ifft(psi)

        # Renormalize
        Norm = sum(abs.(psi).^2) * Transf.dx * Transf.dy * Transf.dz
        psi = sqrt(Params.N) .* psi ./ sqrt(Norm)
        psi = real(psi)
        muchem = chemicalpotential(psi, Params, Transf, VDk, V, HT)

        # Plotting/observables loop
        if mod(t_idx, 1) == 0

            # Energy
            E = energytotal(psi, Params, Transf, VDk, V, Ftherm)
            E = E / Norm
            # push!(Observ.EVec, E)

            # Chemical potential
            # push!(Observ.mucVec, muchem)

            # Normalized residual
            res = norm_resid(psi, Params, Transf, VDk, V, muchem, HT)
            # push!(Observ.residual, res)

            # Dipole-dipole energy
            # Eddi = E_dd(psi, Params, Transf, VDk)
            # push!(Observ.EddiVec, Eddi)

            # Observ.res_idx += 1

            runningplot(psi, Params, Transf)

            # Save plots and params
            # if Params.save_plots == 1
            #     savefig(joinpath(Params.saving_dir, "TeGPE_groundstate.png"), dpi=300)
            # end
            # if Params.save_params == 1
            #     filename = joinpath(Params.saving_dir, "params.txt")
            #     open(filename, "w") do f
            #         println(f, Params)
            #     end
            # end

            # Save data (using JLD2 instead of MATLAB .mat)
            # @save joinpath(Params.saving_dir, "TeGPE_groundstate.jld2") psi muchem Observ t_idx Transf Params VDk V

            # Convergence checks
            # relres = abs(Observ.residual[Observ.res_idx] - Observ.residual[Observ.res_idx - 1]) /
            #          Observ.residual[Observ.res_idx]
            # relEddi = abs((Observ.EddiVec[Observ.res_idx] - Observ.EddiVec[Observ.res_idx - 1]) /
            #               Observ.EddiVec[Observ.res_idx])

            if relEddi < Params.stop_rel_Eddi && Params.stop_rel_Eddi_flag
                println("Relative change in Eddi smaller than $(Params.stop_rel_Eddi), simulation stops.")
                break
            end

            if relres < 1e-4
                if relres < Params.stop_relres && Params.stop_relres_flag == 1
                    println("Relative residual smaller than $(Params.stop_relres), simulation stops.")
                    break
                elseif AdaptIdx > 2 && abs(dt) > Params.mindt
                    dt = dt / 2
                    println("Time step changed to $dt")
                    AdaptIdx = 0
                elseif AdaptIdx > 3 && abs(dt) < Params.mindt
                    break
                else
                    AdaptIdx += 1
                end
            else
                AdaptIdx = 0
            end
        end

        # NaN check
        if any(isnan.(psi))
            println("NaN detected in psi, stopping.")
            break
        end

        t_idx += 1
    end
    h5open("./compdata/ssfm_imag_julia.h5", "w") do file
            write(file, "psi", psi)
            write(file, "V", V)
            write(file, "VDk", VDk)
            write(file, "dx", Transf.dx)
            write(file, "dy", Transf.dy)
            write(file, "dz", Transf.dz)
            write(file, "muchem", muchem)
    end

    # Final observables
    E = energytotal(psi, Params, Transf, VDk, V, Ftherm)
    Norm = sum(abs.(psi).^2) * Transf.dx * Transf.dy * Transf.dz
    E = E / Norm
    push!(Observ.EVec, E)
    push!(Observ.mucVec, muchem)

    res = norm_resid(psi, Params, Transf, VDk, V, muchem, HT)
    push!(Observ.residual, res)

    Observ.res_idx += 1

    return psi, muchem, Observ, t_idx
end