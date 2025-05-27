using Plots
using Plots.Measures

function set_common_plot_style()
    default(
        fontfamily="Computer Modern",
        titlefontsize=19,
        guidefontsize=10,
        tickfontsize=16,
        legendfontsize=16,
        left_margin = 2mm,
        bottom_margin = 2mm,
        right_margin = 2.5mm,
        top_margin = 2.5mm,
        legend=:topright,
        framestyle=:box,
        minorgridalpha=.04,
        gridwidth=0.5,
        gridalpha=0.5,
        labelfontsize=21
    )
end