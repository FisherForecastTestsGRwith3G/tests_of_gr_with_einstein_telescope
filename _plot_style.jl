using Plots
using Plots.Measures
using ColorSchemes

function set_common_plot_style()
    # Default value is scale = 1.0, but e.g. scale = 1.4 is better for side-by-side plots
    scale = 1.0
    default(
        fontfamily="Computer Modern",
        titlefontsize=round(Int, 19 * scale),
        guidefontsize=round(Int, 10 * scale),
        tickfontsize=round(Int, 16 * scale),
        legendfontsize=round(Int, 16 * scale),
        colorbar_titlefontsize = round(Int, 18 * scale),
        labelfontsize=round(Int, 21 * scale),
        left_margin = 8mm,
        right_margin = 2mm,
        bottom_margin = 12mm,
        top_margin = 10mm,
        legend=:topright,
        framestyle=:box,
        minorgridalpha=.04,
        gridwidth=0.5,
        gridalpha=0.5
    )
end

function default_plot_dimensions()
    # Set default plot dimensions
    return (plotHeight=800, plotWidth=1200, plotDpi=300, padding = (8mm, 5mm))
end

function smaller_margins_subplots()
    # Set default plot margins for subplots
    return (left_margin = 2mm, right_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm, left_margin_with_ylabel = 8mm)
end

function labels_from_networks(network_names)
    # Define a dictionary to map network names to their labels
    network_labels = Dict(
        "ETS" => L"\texttt{T}",
        "network_45_15km" => L"\texttt{2L\_45}",
        "network_0_15km" => L"\texttt{2L\_0}",
        "LHV" => L"\texttt{HLV}",
        "LHVK" => L"\texttt{HKLV}",
        "final_run_LVK_LHV" => L"\texttt{HLV\ (PhenomHM)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_LVK_PhenomD_LHV" => L"\texttt{HLV\ (PhenomD)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_ET_ETS" => L"\texttt{T\ (PhenomHM)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_ET_PhenomD_ETS" => L"\texttt{T\ (PhenomD)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_ET_network_45_15km" => L"\texttt{2L\_45\ (PhenomHM)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_ET_PhenomD_network_45_15km" => L"\texttt{2L\_45\ (PhenomD)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_ET_network_0_15km" => L"\texttt{2L\_0\ (PhenomHM)}", # Needed when overloading the detector networks as waveform models in plot C
        "final_run_ET_PhenomD_network_0_15km" => L"\texttt{2L\_0\ (PhenomD)}", # Needed when overloading the detector networks as waveform models in plot C
        "ET_GW150914_like_fmin_2Hz_ETS" => L"GW150914-like, ET $(\texttt{T})$, $f_{\mathrm{min}} = 2$ Hz", # Needed when overloading the detector networks as waveform models in plot C [For GW150914 like event]
        "ET_GW150914_like_fmin_5Hz_ETS" => L"GW150914-like, ET $(\texttt{T})$, $f_{\mathrm{min}} = 5$ Hz", # Needed when overloading the detector networks as waveform models in plot C [For GW150914 like event]
        "ET_GW150914_like_fmin_10Hz_ETS" => L"GW150914-like, ET $(\texttt{T})$, $f_{\mathrm{min}} = 10$ Hz", # Needed when overloading the detector networks as waveform models in plot C [For GW150914 like event]
        "ET_GW150914_like_fmin_15Hz_ETS" => L"GW150914-like, ET $(\texttt{T})$, $f_{\mathrm{min}} = 15$ Hz", # Needed when overloading the detector networks as waveform models in plot C [For GW150914 like event]
        "ET_GW150914_like_fmin_20Hz_ETS" => L"GW150914-like, ET $(\texttt{T})$, $f_{\mathrm{min}} = 20$ Hz" # Needed when overloading the detector networks as waveform models in plot C [For GW150914 like event]
    )

    # Map the network names to their corresponding labels
    return [network_labels[name] for name in network_names]
end

#TODO: Remove this function
function labels_from_PN_orders(PN_orders)
    # Define a dictionary to map PN orders to their labels
    PN_labels = Dict(
        "-1" => L"\varphi_{-2}",
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
    markersize = 8
    
    # Define a color palette
    # palette_var = palette(:seaborn_colorblind)
    palette_var = ColorSchemes.seaborn_colorblind #[:orange, :blue, :green, :purple, :red, :cyan, :magenta, :yellow]    #viridis

    return markers, markersize, palette_var
end

function get_relative_x_offset_network(current_network_index, total_number_networks)
    # Calculate the relative x offset (from -1.0 to 1.0) for each network, used to avoid overlap. I assume current_network_index starts from 1.
    if total_number_networks == 1
        return 0.0
    end
    return ((current_network_index - 1) / (total_number_networks - 1)) * 2.0 - 1.0
end

const labels_from_pn_orders = Dict(
        "-1" => L"\varphi_{-2}",
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

const color_from_pn_orders = Dict(
        "-1" => get_markers_and_palette()[3][1],
        "0" => get_markers_and_palette()[3][2],
        "0.5" => get_markers_and_palette()[3][3],
        "1" => get_markers_and_palette()[3][4],
        "1.5" => get_markers_and_palette()[3][5],
        "2" => get_markers_and_palette()[3][6],
        "log(2.5)" => get_markers_and_palette()[3][7],
        "3" => get_markers_and_palette()[3][8],
        "log(3.)" => get_markers_and_palette()[3][9],
        "3.5" => get_markers_and_palette()[3][10]
    )
