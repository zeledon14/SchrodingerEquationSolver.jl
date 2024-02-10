using SchrodingerEquationSolver
using Test

@testset "SchrodingerEquationSolver.jl" begin
    include("uS1HydrogenTest.jl")
    include("OneDSchrodingerEquationSolverTest.jl")
    include("HarmincOcillatorEigenStatesSearchTest.jl")
    include("uS1HydrogenHartreeTest.jl")
    include("HeTest.jl")
end
