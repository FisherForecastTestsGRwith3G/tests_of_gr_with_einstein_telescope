using HDF5
using TOML
using LaTeXStrings
using CairoMakie

include("_config_parser_grid.jl")
include("../scripts_to_run/_plot_style.jl")
include("../create_single_event_datasets/createSED.jl")
using .createSED: pnoString

# obtain config file from input
config_file_name = ARGS[1]
println("Using config file: $(config_file_name)\n")

configs = read_config_grid(config_file_name)

pn_string = configs["pn_orders"][1]
pn_tag = pnoString(pn_string)
network = configs["network"]

result_folder_name = abspath(joinpath(@__DIR__, configs["outdir"], "grid", network, pn_tag))
plot_folder_name = abspath(joinpath(@__DIR__, configs["plot_outdir"], "grid", network, pn_tag))
mkpath(plot_folder_name)

result_file_name = joinpath(result_folder_name, "results.h5")
isfile(result_file_name) || throw(ArgumentError("Result file not found: $(result_file_name). Run B_grid_plot.jl first."))

mu_vec = configs["mu_vec"]
sigma_vec = configs["sigma_vec"]
res = nothing

h5open(result_file_name, "r") do file
    res = read(file, "results")
    if haskey(file, "mu_vec")
        mu_vec = read(file, "mu_vec")
    end
    if haskey(file, "sigma_vec")
        sigma_vec = read(file, "sigma_vec")
    end
    if haskey(file, "pn_string")
        pn_string = read(file, "pn_string")
    end
end

fig_file_name = joinpath(plot_folder_name, "grid_plot.pdf")
println("\nSaving grid plot in: $(fig_file_name)")

### leave unchanged
fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))
####

tickfontsize = 24
titlefontsize = 26
labelfontsize = 30
lim = collect(0:0.4:3)
lim_label = []
for lim_i in lim
    push!(lim_label, string.(Int.(round.(10 .^ lim_i, sigdigits = 1))))
end

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"\sigma", title = createSED.pnoLatex(pn_string))

# Create a heatmap with the results
hm = Makie.heatmap!(ax, mu_vec, sigma_vec, log10.(res), colormap = FIG9_IMPROVEMENT_COLORMAP)
Makie.Colorbar(fig[1, 2], hm, size = 20,
    ticklabelsize = 20, ticks = (lim, lim_label))

Makie.Label(fig[1, 2, Top()], L"n_{\text{events}}", fontsize = labelfontsize)

ax.titlesize = titlefontsize
ax.xlabelsize = labelfontsize
ax.ylabelsize = labelfontsize

ax.xticklabelsize = tickfontsize
ax.yticklabelsize = tickfontsize
ax.yticklabelspace = 50.

Makie.save(fig_file_name, fig)
