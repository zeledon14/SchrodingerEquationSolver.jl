using SchrodingerEquationSolver
using Test

@testset "SchrodingerEquationSolver.jl" begin
    include("CTest.jl")
    include("Eigenstates_search_1-D_gaussian_potentialTest.jl")
    include("Eigenstates_search_1-D_quantum_harmonic_oscillator_linear_gridTest.jl")
end

#to run test

#= in SchrodingerEquationSolver directory

julia

using Pkg
Pkg.activate(".")
Pkg.test() =#

