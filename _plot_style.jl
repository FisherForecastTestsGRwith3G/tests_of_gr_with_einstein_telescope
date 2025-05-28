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

function default_plot_dimensions()
    # Set default plot dimensions
    return (plotHeight=500, plotHeiplotWidthght=1200, plotDpi=300)
end

function labels_from_networks(network_names)
    # Define a dictionary to map network names to their labels
    network_labels = Dict(
        "ETS" => L"\textbf{T}",
        "network_45_15km" => L"\textbf{2L\_45}",
        "network_0_15km" => L"\textbf{2L\_0}",
    )

    # Map the network names to their corresponding labels
    return [network_labels[name] for name in network_names]
end

function labels_from_PN_orders(PN_orders)
    # Define a dictionary to map PN orders to their labels
    PN_labels = Dict(
        "-1" => L"\varphi_{-1}",
        "0" => L"\varphi_{0}",
        "0.5" => L"\varphi_{1}",
        "1" => L"\varphi_{2}",
        "1.5" => L"\varphi_{3}",
        "2" => L"\varphi_{4}",
        "log(2.5)" => L"\varphi_{5\,\ell}",
        "3" => L"\varphi_{6}",
        "log(3.)" => L"\varphi_{6\,\ell}",
        "3.5" => L"\varphi_{7}"
    )

    # Map the PN orders to their corresponding labels
    return [PN_labels[order] for order in PN_orders]
end

function get_markers_and_palette()
    # Define a list of different markers to enhance readability
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :hexagon]
    
    # Define a color palette
    palette = palette(:seaborn_colorblind) #[:orange, :blue, :green, :purple, :red, :cyan, :magenta, :yellow]    #viridis

    return markers, palette
end