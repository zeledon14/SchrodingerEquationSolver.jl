using AutomaticDocstrings
using Test

import .SchrodingerEquationSolver as ses
import .ses.Grids as Grids
import .ses.Potentials as Potentials
import .ses.MathUtils as MathUtils
import .ses.Hydrogen as Hydrogen
import .ses.InitialConditions as InitialConditions
import .ses.OneDSchrodingerEquationSolver as OneDSchrodingerEquationSolver
import .ses.OneDPoissonEquationSolver as OneDPoissonEquationSolver


import .ses.EigenvalueFinders as EigenvalueFinders
import .ses.AtomBasisSet as AtomBasisSet
import .ses.Density as Density
import .ses.ExchangeCorrelation as ExchangeCorrelation
@testset "HeTest" begin
        
    r_max::Float64=50.0;
    Z::Int64=2;
    grid= ses.Grids.exponential_grid(r_max, Z, 0.001504);
    grid_sqrt= grid.^2.0;
    N= size(grid)[1];

    #init potentials
    V_colu::Vector{Float64}= Potentials.coulomb_potential(Z, grid);
    V_hartree::Vector{Float64}=zeros(Float64, N);
    V_x::Vector{Float64}=zeros(Float64, N);
    V_c::Vector{Float64}=zeros(Float64, N);
    E_xp::Vector{Float64}=zeros(Float64, N);
    E_cp::Vector{Float64}=zeros(Float64, N);
    V_xcp::Vector{Float64}=zeros(Float64, N);
    density_in::Vector{Float64}= zeros(Float64, N);

    E_total::Float64=1.0
    E_total_before::Float64=2.0;
    E_eigen::Float64=0.0
    V_xc::Float64= 0.0
    E_hartree::Float64= 0.0
    E_x::Float64= 0.0
    E_c::Float64= 0.0
#    C_in= 0.0
#    C_out= 0.0
    basis= AtomBasisSet.init_atom_basis_set(Z, grid);

    while abs(E_total - E_total_before) > 10.0e-12
        E_eigen=0.0
        for i_orbi in basis.orbitals
            V_angu= Potentials.angular_potential(i_orbi.l, grid);
            V_effe= V_colu .+ V_angu .+ V_hartree .+ V_x .+ V_c;
    
            E_grid= Grids.uniform_grid(i_orbi.E - 0.5*i_orbi.E, i_orbi.E + 0.5*i_orbi.E, 150);
    
            E_intervals, bad_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, V_effe,
                                                            grid,InitialConditions.atom, i_orbi.l);

            u_temp, ei_temp= EigenvalueFinders.illinois_eigenvalue_finder(E_intervals[1], V_effe, grid,InitialConditions.atom);
            i_orbi.E=ei_temp
            i_orbi.u=u_temp
            E_eigen+= i_orbi.occu*ei_temp
    
        end
    
        E_total_before= float(E_total)
        density_out= Density.calculate_density(basis);
        density_in= Density.linear_mixing(density_in, density_out, alpha=0.25);
        V_hartree= OneDPoissonEquationSolver

.solver(Z, density_in, grid);
        V_x, E_xp, V_c, E_cp= ExchangeCorrelation.potentials(density_in);
    
        V_xcp= V_x .+ V_c;
        V_xc= (4.0*pi)*MathUtils.integral((V_xcp.*density_in.*grid_sqrt), grid)
        E_hartree= 0.5*(4.0*pi)*(MathUtils.integral((V_hartree.*density_in.*grid_sqrt), grid))
        E_x= (4.0*pi)*MathUtils.integral((E_xp.*density_in.*grid_sqrt), grid)
        E_c= (4.0*pi)*MathUtils.integral((E_cp.*density_in.*grid_sqrt), grid)
        #C_in= (4.0*pi)*MathUtils.integral((density_in.*grid_sqrt), grid)
        #C_out= (4.0*pi)*MathUtils.integral((density_out.*grid_sqrt), grid)
        E_total= E_eigen - E_hartree + E_x + E_c - V_xc
    end
    @test abs(E_total - -2.834836) < 5.0e-7;
end