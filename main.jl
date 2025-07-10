#main.jl

using Pkg
Pkg.activate(@__DIR__)

using SchrodingerEquationSolver
using SchrodingerEquationSolver: libxc_DFTAtom, AtomBasisSet, ExchangeCorrelation
using Plots
using CSV
using DataFrames
using PrettyTables
using Libxc
using YAML

function main(config_path::String)
    println("Loading config from: $config_path")
    config = YAML.load_file(config_path)


    elements_list = config["elements"]
    r_max_list = config["r_max"]
    if length(elements_list) > 0 && length(elements_list) == length(r_max_list)
        println("Everyone of the $(length(elements_list)) elements has its own r_max")
    elseif length(elements_list) > length(r_max_list) && length(r_max_list) == 1
        r_max_list = fill(r_max_list[1], length(elements_list))
        println("Using the same r_max for all the $(length(elements_list)) elements")
    else
        error("The r_max or the elements list is an empty or mismatched array")
    end
    potential_type = config["potential_type"]
    s= get(config, "s_blum_potential", 200.0);
    r_onset= get(config, "r_onset", 4.0);


    println("Elements: $elements_list")
    println("r_max: $r_max_list")
    println("Potential type: $potential_type")

    for (r_max_i, z_i) in zip(r_max_list, elements_list)
        out = libxc_DFTAtom.calculate_atomic_basis_set(z_i, r_max=r_max_i, s=s, r_onset= r_onset);
    end

    # You can now call your basis set function here
end

# Parse command-line arguments
function parse_args()
    for (i, arg) in enumerate(ARGS)
        if arg == "--config_path" && i < length(ARGS)
            return ARGS[i+1]
        end
    end
    error("Missing --config_path argument.")
end

# Entry point
if abspath(PROGRAM_FILE) == @__FILE__
    config_path = parse_args()
    main(config_path)
end