# Shared color palette
const COLORS = [:royalblue, :crimson, :forestgreen, :darkorange, :purple, :teal]

"""
    save_figure(fig, path_stem)

Save figure as both PNG and PDF. `path_stem` should not include extension.
Example: `save_figure(fig, "figures/sigma_z_comparison")`
"""
function save_figure(fig, path_stem::String)
    mkpath(dirname(path_stem))
    CairoMakie.save(path_stem * ".png", fig; px_per_unit=3)
    CairoMakie.save(path_stem * ".pdf", fig)
end

"""
    plot_dynamics(results::Vector{SimulationResult}; observable="sigma_z", title="")

Plot one observable vs time for multiple SimulationResults on the same axes.
"""
function plot_dynamics(results::Vector{SimulationResult};
                       observable::String="sigma_z",
                       title::String="")
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1];
        xlabel="t",
        ylabel="⟨$(observable)⟩",
        title=title
    )
    for (i, r) in enumerate(results)
        lines!(ax, r.t, r.observables[observable];
            color=COLORS[mod1(i, length(COLORS))],
            label=r.method
        )
    end
    axislegend(ax; position=:rt)
    return fig
end

"""
    plot_comparison(r1::SimulationResult, r2::SimulationResult; observable="sigma_z")

Plot two results overlaid plus their difference in a subplot below.
"""
function plot_comparison(r1::SimulationResult, r2::SimulationResult;
                         observable::String="sigma_z")
    fig = Figure(size=(600, 600))

    # Top: overlaid dynamics
    ax1 = Axis(fig[1, 1]; ylabel="⟨$(observable)⟩", title="$(r1.method) vs $(r2.method)")
    lines!(ax1, r1.t, r1.observables[observable]; color=COLORS[1], label=r1.method)
    lines!(ax1, r2.t, r2.observables[observable]; color=COLORS[2], label=r2.method, linestyle=:dash)
    axislegend(ax1; position=:rt)

    # Bottom: difference
    t_common = r1.t
    diff = r1.observables[observable] .- r2.observables[observable]
    ax2 = Axis(fig[2, 1]; xlabel="t", ylabel="Δ⟨$(observable)⟩")
    lines!(ax2, t_common, diff; color=:black)
    hlines!(ax2, [0.0]; color=:gray, linestyle=:dash)

    return fig
end

"""
    plot_convergence(results::Vector{SimulationResult}, param_name::String, param_values;
                     observable="sigma_z")

Plot convergence: multiple curves for different values of a convergence parameter.
"""
function plot_convergence(results::Vector{SimulationResult},
                          param_name::String,
                          param_values;
                          observable::String="sigma_z")
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1];
        xlabel="t",
        ylabel="⟨$(observable)⟩",
        title="Convergence: $(param_name)"
    )
    for (i, (r, v)) in enumerate(zip(results, param_values))
        lines!(ax, r.t, r.observables[observable];
            color=COLORS[mod1(i, length(COLORS))],
            label="$(param_name)=$(v)"
        )
    end
    axislegend(ax; position=:rt)
    return fig
end
