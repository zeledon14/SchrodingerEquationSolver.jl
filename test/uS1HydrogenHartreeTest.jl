using AutomaticDocstrings
using Test

import .SchrodingerEquationSolver as ses
import .ses.Grids as Grids
import .ses.Potentials as Potentials
import .ses.MathUtils as MathUtils
import .ses.Hydrogen as Hydrogen
import .ses.InitialConditions as InitialConditions
import .ses.OneDSchrodingerEquationSolver as odses
import .ses.OneDPoissonEquationSolver as odpes
import .ses.EigenvalueFinders as EigenvalueFinders
import .ses.AtomBasisSet as AtomBasisSet
import .ses.Density as Density

@testset "uS1HydrogenHartreeTest" begin
        
    r_max::Float64=10.0;
    Z::Int64=1;
    grid= ses.Grids.exponential_grid(r_max, Z);
    N::Int64=size(grid)[1];
    basis= AtomBasisSet.init_atom_basis_set(Z, grid);
    v_colu= Potentials.coulomb_potential(Z, grid);
    for i_orbi in basis.orbitals
        v_angu= Potentials.angular_potential(i_orbi.l, grid);
        v_effe= v_colu + v_angu;

        E_grid= Grids.uniform_grid(i_orbi.E - 0.1*i_orbi.E, i_orbi.E + 0.05*i_orbi.E, 50);

        E_intervals, bad_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, v_effe,
                                                        grid,InitialConditions.atom, i_orbi.l);
        u_temp, ei_temp= EigenvalueFinders.illinois_eigenvalue_finder(E_intervals[1], v_effe, grid,InitialConditions.atom);
        i_orbi.E=ei_temp
        i_orbi.u=u_temp

    end
    density= Density.calculate_density(basis);
    V_hartree= odpes.solver(Z, density, grid);
    U_hartree= V_hartree.*grid;
    U_hartree_target=Hydrogen.U_hartree(grid);
    error= MathUtils.error_difference(U_hartree,U_hartree_target)
    @test error < 1.000e-7;
end