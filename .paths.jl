## Not so sure how you want to implement the different path for each of us, 
# so that there is no need to change the code for each of us.

function whoIsThere(PhD)
    path_catalog = ""
    path_output = ""
    if PhD == "Andrea"
        println("Andrea is here")
        # my paths are
        path_catalog = "../gwbeast/catalogs"
        path_output = "../gwbeast/"
    elseif PhD == "Matteo"
        println("Matteo is here")
        # my paths are

    elseif PhD == "Joachim"
        println("Joachim is here")
        # my paths are
    else
        println("I don't know who you are")

    end

    return path_catalog, path_output
end
