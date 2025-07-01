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

function main()
    arra = [2, 6, 8, 10, 14, 18, 36, 54, 86]
    r_max = 50.0

    out = [libxc_DFTAtom.calculate_atomic_basis_set(Z, r_max=r_max) for Z in arra]

    # Show basic info from output (optional)
    for (Z, basis) in zip(arra, out)
        println("Z = $Z, Number of radial functions = ", length(basis.radial_functions))
    end

    # You could also optionally save the data or plot it
end

main()