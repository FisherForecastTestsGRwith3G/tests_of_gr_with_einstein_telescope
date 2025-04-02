import JSON

const pn_order_dic = Dict(
    "-1"       => (-1.0    ,"minus_one"),
    "0"        => (0.0     ,"zero"), 
    "0.5"      => (0.5     ,"half"),
    "1"        => (1.0     ,"one"), 
    "1.5"      => (1.5     ,"one_half"),
    "2"        => (2.0     ,"two"),
    "log(2.5)" => (log(2.5),"log_two_half"),
    "3"        => (3.0     ,"three"), 
    "log(3.)"  => (log(3.) ,"log_three"),
    "3.5"      => (3.5     ,"three_half")
);

function getUserConfigs()

    user_configs = open("config_files/user_configs.json","r") do f
        user_configs = JSON.parse(f)
    end

    user = user_configs["user"]
    println("Logged in as user: $(user)")

    return user_configs
end

function readConfigForA(json_file_name)

    config_dic = open(json_file_name,"r") do f
        config_dic = JSON.parse(f)
    end
    
    configs = merge(config_dic["global"], config_dic["a_specific"])
    run_tag = config_dic["header"]
    return configs, run_tag
end

function readConfigForB(json_file_name)

    config_dic = open(json_file_name,"r") do f
        config_dic = JSON.parse(f)
    end
    configs = merge(config_dic["global"], config_dic["b_specific"])
    configs = merge(configs, config_dic["mcmc_settings"])

    # set defaul values for the limits
    for key in ["mu_lims", "sigma_lims"]
        for pno in keys(pn_order_dic)
            
            tmp = configs[key][pno]
            if tmp == []
                configs[key][pno] = nothing
            else
                configs[key][pno] = (Float64(tmp[1]), Float64(tmp[2]))
            end
        end
    end

    # add injected values from injected mu and sigma
    # TODO: Need to change the way the injection works in the future
    #       Then, in the future we will have to also adopt this 
    configs["injected_values"] = Dict{String, Tuple{Float64,Float64}}()
    for pno in keys(pn_order_dic)
        configs["injected_values"][pno] = (
            Float64(config_dic["a_specific"]["mu"]),
            Float64(config_dic["a_specific"]["sigma"]),
            )
    end

    run_tag = config_dic["header"]

    return configs, run_tag
end
