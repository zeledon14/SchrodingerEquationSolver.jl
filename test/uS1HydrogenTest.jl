using AutomaticDocstrings
using Test

import .SchrodingerEquationSolver.Grids as Grids
import .SchrodingerEquationSolver.Potentials as Potentials
import .SchrodingerEquationSolver.Hydrogen as Hydrogen
import .SchrodingerEquationSolver.IntegralNumericalMethods as IntegralNumericalMethods
import .SchrodingerEquationSolver.MathUtils as MathUtils
@testset "uS1HydrogenTest" begin
        
        r_max::Float32=10.0;
        l::Int32=0;
        Z::Int32=1;
        E::Float32= -0.5000;

        grid= Grids.exponential_grid(r_max, Z);
        v_colu= Potentials.coulomb_potential(Z, grid);
        u_s1_hydr= Hydrogen.u_s1_hydrogen(grid);
        init_valu1_fwrd::Float32=u_s1_hydr[1];#grid[1]^(l+1.0)
        init_valu2_fwrd::Float32=u_s1_hydr[2];#grid[2]^(l+1.0)
        init_valu1_bwrd::Float32=u_s1_hydr[end];#grid[1]^(l+1.0)
        init_valu2_bwrd::Float32=u_s1_hydr[end-1];#grid[2]^(l+1.0)
        v_effe= v_colu;
        f::Vector{Float32}= 2.0.*(v_effe .- E);
        #given f find the classical turning points
        g=zeros(Float32, size(f)[1]);
        turn_pnts= MathUtils.turning_points_indices(f);

        #do forward integration of radial shcrodinger equation u
        u_fwd= IntegralNumericalMethods.integrate_second_order_DE(grid,g,f,
                                                                    init_valu1_fwrd,init_valu2_fwrd);
        #do backward integreation of the radial shcrodinger equation u 
        u_bwd= reverse(IntegralNumericalMethods.integrate_second_order_DE(reverse(grid),g,reverse(f),
                        init_valu1_bwrd,init_valu2_bwrd));
        #rescale u_fwd, u_bwd to make u_fwd[turn_pnts[1]] = u_bwd[turn_pnts[1]]
        u_fwd, u_bwd= MathUtils.rescale!(u_fwd, u_bwd, turn_pnts[1]);
        #merge solutions
        u_merged, merge_value= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
        temp= ((u_merged .- u_s1_hydr).^2).^0.5;
        @test sum(temp)/length(temp) < 1.325e-7;
end