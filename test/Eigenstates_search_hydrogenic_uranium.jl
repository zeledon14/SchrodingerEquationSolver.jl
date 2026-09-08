using AutomaticDocstrings
using Test

using SchrodingerEquationSolver
using SchrodingerEquationSolver:   Potentials, MathUtils, Hydrogen, InitialConditions,
                                   IntegralNumericalMethods, OneDSchrodingerEquationSolver,
                                   EigenvalueFinders, AtomicBasisSets
@testset "Eigenstates_search_hydrogenic_uranium" begin

    Z=92;
    r_max= 50.0;
    grid_stru= Grids.ExponentialGrid(r_max, Z);
    grid= grid_stru.grid;
    r_min=grid_stru.grid[1];
    r_max=grid_stru.grid[end];
    grid_i=grid_stru.grid_i;
    uran_basis= AtomicBasisSets.AtomBasisSet(Z, grid);
    l_max= maximum([elem.l for elem in uran_basis.orbitals])
    lns=Dict{Int64,Dict}(l=>Dict{String,Int64}("n_max"=>1, "n_min"=>l+1) for l in (0:l_max))

    initial_condition_function_x_min = InitialConditions.atom_like_poly_at_r_min;
    initial_condition_function_x_max = InitialConditions.exponential_decay_at_r_ref;
    for (i, iorbi) in enumerate(uran_basis.orbitals)
        l= iorbi.l;
        n=iorbi.n;
        E_target=iorbi.E
        v_effe =(
            Potentials.angular_potential(l, grid) .+
            Potentials.coulomb_potential(Z, grid));
        
        E_pred_inte= (1.03*E_target, 0.973*E_target);

        u,E_pred= EigenvalueFinders.illinois_eigenvalue_finder(E_pred_inte,
        v_effe, grid_stru, 
        initial_condition_function_x_min,
        initial_condition_function_x_max,l)
        @test abs(E_pred - E_target) < 10.0e-9;
    end
    
end