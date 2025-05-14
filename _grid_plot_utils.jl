


# function modify_configs(config_name, header, mu, sigma)
    
#     open(config_name, "w") do

#     # Modify the configs dictionary to set the values of mu and sigma
#     configs["mu"] = mu
#     configs["sigma"] = sigma
#     configs["header"] = header
#     # configs["simulation_tag"] = header * "_mu_" * string(mu) * "_sigma_" * string(sigma)

#     return configs
# end


# using JSON
function modify_configs(config_file_name, config_file_name_out, header, mu, sigma, PN, network_list, n_events)

    # Read the JSON file into a Julia dictionary
    config_dic = open(config_file_name,"r") do f
        config_dic = JSON.parse(f)
    end

    config_glob =config_dic["global"]
    config_a = config_dic["a_specific"]

    # Modify the values
    config_a["mu"] = mu
    config_a["sigma"] = sigma
    config_a["n_events"] = n_events

    config_dic["header"] = header
    config_glob["pn_waveforms"] = [PN]
    config_glob["network_list"] = network_list


    # Write the modified dictionary back to the file
    open(config_file_name_out, "w") do io
        JSON.print(io, config_dic, 4)
    end
end