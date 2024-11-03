using AutomaticDocstrings
using Test

import .SchrodingerEquationSolver.Grids as Grids
import .SchrodingerEquationSolver.Potentials as Potentials
import .SchrodingerEquationSolver.Hydrogen as Hydrogen
import .SchrodingerEquationSolver.IntegralNumericalMethods as IntegralNumericalMethods
import .SchrodingerEquationSolver.MathUtils as MathUtils
import .SchrodingerEquationSolver.InitialConditions as InitialConditions
import .SchrodingerEquationSolver.OneDSchrodingerEquationSolver as OneDSchrodingerEquationSolver
@testset "OneDSchrodingerEquationSolverTest" begin
        
        r_max::Float64=10.0;
        Z::Int64=1;
        l::Int64=0;
        E::Float64= -0.5000;

        grid= Grids.exponential_grid(r_max, Z);
        v_colu= Potentials.coulomb_potential(Z, grid);
        u_s1_hydr= Hydrogen.u_s1_hydrogen(grid);
        u_s1_hydr_norm= MathUtils.normalize!(u_s1_hydr, grid);
        init_valu1_fwrd, init_valu2_fwrd,
        init_valu1_bwrd, init_valu2_bwrd =InitialConditions.atom(grid, E, l);
        v_effe= v_colu;
        u_merged, merge_value= OneDSchrodingerEquationSolver.solver(E,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid);
        error= MathUtils.error_difference(u_merged,u_s1_hydr_norm)
        @test error < 1.60e-7
end