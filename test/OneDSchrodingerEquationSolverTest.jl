using AutomaticDocstrings
using Test

import .SchrodingerEquationSolver.Grids as Grids
import .SchrodingerEquationSolver.Potentials as Potentials
import .SchrodingerEquationSolver.Hydrogen as Hydrogen
import .SchrodingerEquationSolver.IntegralNumericalMethods as IntegralNumericalMethods
import .SchrodingerEquationSolver.MathUtils as MathUtils
import .SchrodingerEquationSolver.OneDSchrodingerEquationSolver as odses
@testset "OneDSchrodingerEquationSolverTest" begin
        
        r_max::Float32=10.0;
        Z::Int32=1;
        E::Float32= -0.5000;

        grid= Grids.exponential_grid(r_max, Z);
        v_colu= Potentials.coulomb_potential(Z, grid);
        u_s1_hydr= Hydrogen.u_s1_hydrogen(grid);
        u_s1_hydr_norm= MathUtils.normalize!(u_s1_hydr, grid);
        init_valu1_fwrd::Float32=u_s1_hydr[1];#grid[1]^(l+1.0)
        init_valu2_fwrd::Float32=u_s1_hydr[2];#grid[2]^(l+1.0)
        init_valu1_bwrd::Float32=u_s1_hydr[end];#grid[1]^(l+1.0)
        init_valu2_bwrd::Float32=u_s1_hydr[end-1];#grid[2]^(l+1.0)
        v_effe= v_colu;
        u_merged, merge_value= odses.solver(E,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid);
        temp= ((u_merged .- u_s1_hydr_norm).^2).^0.5;
        @test sum(temp)/length(temp) < 1.39
end