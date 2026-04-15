using CairoMakie
using DelimitedFiles
using LaTeXStrings

const ASD_DATA_DIR = joinpath(@__DIR__, "..", "create_single_event_datasets", "psd_data", "hlv_curves")
const OUTPUT_DIR = joinpath(@__DIR__, "results")

const DETECTOR_COLORS = Dict(
    :L1 => "#4ba6ff",
    :H1 => "#ee0000",
    :V1 => "#9b59b6",
)

const ASD_FILES = Dict(
    "O3a" => Dict(
        :L1 => joinpath(ASD_DATA_DIR, "O3a", "O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_strain_asd.txt"),
        :H1 => joinpath(ASD_DATA_DIR, "O3a", "O3-H1-C01_CLEAN_SUB60HZ-1251752040.0_sensitivity_strain_asd.txt"),
        :V1 => joinpath(ASD_DATA_DIR, "O3a", "O3-V1_sensitivity_strain_asd.txt"),
    ),
    "O3b" => Dict(
        :L1 => joinpath(ASD_DATA_DIR, "O3b", "O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt"),
        :H1 => joinpath(ASD_DATA_DIR, "O3b", "O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt"),
        :V1 => joinpath(ASD_DATA_DIR, "O3b", "O3-V1-1265246178_sensitivity_strain_asd.txt"),
    ),
    "O4" => Dict(
        :L1 => joinpath(ASD_DATA_DIR, "O4_preO4estimates", "aligo_O4high.txt"),
        :H1 => joinpath(ASD_DATA_DIR, "O4_preO4estimates", "aligo_O4high.txt"),
        :V1 => joinpath(ASD_DATA_DIR, "O4_preO4estimates", "avirgo_O4high_NEW.txt"),
    ),
)

function load_asd_curves()
    return Dict(
        run => Dict(detector => Matrix{Float64}(loadtxt(path)) for (detector, path) in files)
        for (run, files) in ASD_FILES
    )
end

function loadtxt(path::AbstractString)
    return readdlm(path, Float64)
end

function comparison_title(run_a::String, run_b::String)
    if run_a == "O4"
        return "O4 ASD estimate vs $(run_b)"
    elseif run_b == "O4"
        return "O4 ASD estimate vs $(run_a)"
    end
    return "$(run_a) ASD vs $(run_b)"
end

function run_maturity(run::String)
    return Dict("O3a" => 1, "O3b" => 2, "O4" => 3)[run]
end

function run_alpha(run::String, run_a::String, run_b::String)
    less_mature_run = run_maturity(run_a) < run_maturity(run_b) ? run_a : run_b
    return run == less_mature_run ? 0.5 : 0.95
end

function add_comparison_panel!(fig, grid_position, curves, run_a::String, run_b::String)
    ax = Axis(
        fig[grid_position...];
        xscale=log10,
        yscale=log10,
        xlabel=L"f\,[\mathrm{Hz}]",
        ylabel=L"\mathrm{ASD}\,[1/\sqrt{\mathrm{Hz}}]",
        title=comparison_title(run_a, run_b),
        titlesize=24,
        xlabelsize=22,
        ylabelsize=22,
        xticklabelsize=16,
        yticklabelsize=16,
        xgridvisible=true,
        ygridvisible=true,
    )

    for run in (run_a, run_b)
        for detector in (:L1, :H1, :V1)
            data = curves[run][detector]
            lines!(
                ax,
                view(data, :, 1),
                view(data, :, 2);
                color=(DETECTOR_COLORS[detector], run_alpha(run, run_a, run_b)),
                linestyle=:solid,
                linewidth=3,
                label="$(detector) $(run)",
            )
        end
    end

    xlims!(ax, 10, 4_000)
    ylims!(ax, 1e-24, 1e-18)
    axislegend(ax; position=:rb, framevisible=false, labelsize=14, rowgap=2)
    return ax
end

function main()
    mkpath(OUTPUT_DIR)
    curves = load_asd_curves()

    fig = Figure(size=(2200, 700))
    comparisons = [("O4", "O3a"), ("O4", "O3b"), ("O3a", "O3b")]

    for (idx, (run_a, run_b)) in enumerate(comparisons)
        add_comparison_panel!(fig, (1, idx), curves, run_a, run_b)
    end

    output_base = joinpath(OUTPUT_DIR, "hlv_asd_comparisons")
    save(output_base * ".png", fig)
    save(output_base * ".pdf", fig)
    return nothing
end

main()
