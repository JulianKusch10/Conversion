using Plots

function runningplot(psi, Params, Transf)

    # Density
    n = abs.(psi).^2

    nxz = dropdims(sum(n .* Transf.dy, dims=2), dims=2)
    nyz = dropdims(sum(n .* Transf.dx, dims=1), dims=1)
    nxy = dropdims(sum(n .* Transf.dz, dims=3), dims=3)

    p1 = heatmap(Transf.x, Transf.z, nxz',
        xlabel=L"x [\mu m]", ylabel=L"z [\mu m]", colorbar=false)

    p2 = heatmap(Transf.y, Transf.z, nyz',
        xlabel=L"y [\mu m]", ylabel=L"z [\mu m]", colorbar=false)

    p3 = heatmap(Transf.x, Transf.y, nxy',
        xlabel=L"x [\mu m]", ylabel=L"y [\mu m]", colorbar=false)

    p4 = plot(-log10.(Observ.residual),
        xlabel="steps", ylabel=L"-\log_{10}(r)", legend=false)

    if Params.stop_rel_Eddi_flag == 1
        p5 = plot(Observ.EVec, label="E", xlabel="steps", ylabel=L"E")
        plot!(p5, Observ.EddiVec, label="E_dd", ylabel=L"E_{dd}")
    else
        p5 = plot(Observ.EVec,
            xlabel="steps", ylabel=L"E", legend=false)
    end

    p6 = plot(Observ.mucVec,
        xlabel="steps", ylabel=L"\mu", legend=false)

    title_str = "N = $(Params.N), T = $(Params.T) nK, a = $(Params.as_a0)a₀, θ = $(Params.theta_deg)°"

    plot(p1, p2, p3, p4, p5, p6,
        layout=(2,3),
        size=(1000,600),
        plot_title=title_str)
end