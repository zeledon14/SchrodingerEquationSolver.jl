module SchrodingerEquationSolver

include("Grids.jl");
include("MathUtils.jl");
include("Potentials.jl");
include("IntegralNumericalMethods.jl");
include("OneDSchrodingerEquationSolver.jl");
include("OneDPoissonEquationSolver.jl");
include("InitialConditions.jl");
include("EigenvalueFinders.jl");
include("AtomBasisSet.jl");
include("Density.jl");
include("ExchangeCorrelation.jl")
include("Hydrogen.jl");
end
