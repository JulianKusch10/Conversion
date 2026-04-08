using Plots
using LaTeXStrings
using Printf
using Measures
default(
    guidefontfamily = "Computer Modern",
    tickfontfamily = "Computer Modern",
    legendfontfamily = "Computer Modern",
    titlefontfamily = "Computer Modern"
)
function runningplot(psi, Params, Transf, Observ, save_index)

    # Density
    n = abs.(psi).^2

    nxz = dropdims(sum(n .* Transf.dy, dims=2), dims=2)
    nyz = dropdims(sum(n .* Transf.dx, dims=1), dims=1)
    nxy = dropdims(sum(n .* Transf.dz, dims=3), dims=3)
    N_steps = length(Observ.residual)
    steps = 1:1:N_steps
    steps *= save_index

    p1 = heatmap(Transf.x, Transf.z, nxz',
        xlabel=L"x [\mu m]", ylabel=L"z [\mu m]", colorbar=false)

    p2 = heatmap(Transf.y, Transf.z, nyz',
        xlabel=L"y [\mu m]", ylabel=L"z [\mu m]", colorbar=false)

    p3 = heatmap(Transf.x, Transf.y, nxy',
        xlabel=L"x [\mu m]", ylabel=L"y [\mu m]", colorbar=false)

    p4 = plot(steps, -log10.(Observ.residual),
        xlabel="steps", ylabel=L"-\log_{10}(r)", legend=false, xlims = (0, length(Observ.residual)*save_index))

    p5 = plot(steps, Observ.EVec,
        xlabel="steps", ylabel=L"E", legend=false, xlims = (0, length(Observ.EVec)*save_index))

    p6 = plot(steps, Observ.mucVec,
        xlabel="steps", ylabel=L"\mu", legend=false, xlims = (0, length(Observ.mucVec)*save_index))

    title_str = latexstring(
                "N = ", @sprintf("%.0f", Params.N),
                ",\\; T = ", @sprintf("%.0f", Params.T), "\\,\\mathrm{nK},\\; a = ",
                @sprintf("%.1f", Params.as_a0), "a_0,\\; \\theta = ",
                @sprintf("%.1f", Params.theta_deg), "^\\circ"
)

    p = plot(p1, p2, p3, p4, p5, p6,
        layout=(2,3),
        size=(1000,600),
        plot_title=title_str, 
        left_margin = 5mm, 
        right_margin = 5mm,
        top_margin = 5mm,
        bottom_margin = 5mm)
    savefig(p, "groundstate.png")
end