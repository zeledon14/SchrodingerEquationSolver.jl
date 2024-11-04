using AutomaticDocstrings
using Test

using SchrodingerEquationSolver
using SchrodingerEquationSolver: Grids, Potentials, MathUtils, Hydrogen, InitialConditions,
                                 OneDSchrodingerEquationSolver, OneDPoissonEquationSolver,
                                 EigenvalueFinders, AtomBasisSet, Density, ExchangeCorrelation
@testset "CTest" begin
        
    #Define parameters and produce an exponential grid.
    r_max::Float64=50.0;#Max radius of space grid.
    Z::Int64=6;#Atomic number, also used as the charge of coulomb potential.

        #grid definition
    grid_stru= Grids.init_exponential_grid_structure(r_max, Z);
    N=grid_stru.N;
    @test N == 7835;

    #init potentials
    #Initialization of potentials and energies

    #Initializing coulomb potential due to nuclei charge.
    V_colu::Vector{Float64}= Potentials.coulomb_potential(Z, grid_stru.grid);
    #Initializing Hartree potential due to electron density. 
    V_hartree::Vector{Float64}=zeros(Float64, N);
    #Initializing exchange potential.
    V_x::Vector{Float64}=zeros(Float64, N);
    #Initializing correlation potential.
    V_c::Vector{Float64}=zeros(Float64, N);
    #Initializing energy exchange potential.
    E_xp::Vector{Float64}=zeros(Float64, N);
    #Initializing energy correlation potential.
    E_cp::Vector{Float64}=zeros(Float64, N);
    #Initializing exchange + correlation potential.
    V_xcp::Vector{Float64}=zeros(Float64, N);
    #Initializing the density.
    density_in::Vector{Float64}= zeros(Float64, N);
    #Initializing total energy
    E_total::Float64=1.0;
    #Initializing total energy step before
    E_total_before::Float64=2.0;
    #Initializing energy from energy eigenvalues
    E_eigen::Float64=0.0;
    #Initializing exchange correlation potential value after integral
    V_xc::Float64= 0.0;
    #Initializing Hartree energy
    E_hartree::Float64= 0.0;
    #Initializing exchange energy
    E_x::Float64= 0.0;
    #Initializing correlation enrgy
    E_c::Float64= 0.0;
    #Initializing basis set data structure
    basis= AtomBasisSet.init_atom_basis_set(Z, grid_stru.grid);

    #Energy minimization loop 
    while abs(E_total - E_total_before) > 10.0e-8
        E_eigen=0.0;
        #Loop over every orbital to solve independent particle Schrodinger equation.
        for i_orbi in basis.orbitals
            #println(i_orbi.name)
            #angular potential for l orbital
            V_angu= Potentials.angular_potential(i_orbi.l, grid_stru.grid);
            #Assemble effective potential.
            V_effe= V_colu .+ V_angu .+ V_hartree .+ V_x .+ V_c;
            V_effe_max= maximum(V_effe)
            V_effe_min= minimum(V_effe)
            energy_interval= EigenvalueFinders.guess_energy_interval(i_orbi.E, V_effe_max, V_effe_min);
            E_grid= Grids.uniform_grid(energy_interval[1], energy_interval[2], 200); #List with the energy grid points.

            E_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, V_effe, grid_stru,
                        InitialConditions.atom_exponential_grid,
                            OneDSchrodingerEquationSolver.solver_exponential_grid, numb_inter=1);
            #print("here ")
            u_temp, ei_temp= EigenvalueFinders.illinois_eigenvalue_finder(E_intervals[1], V_effe, 
            grid_stru,InitialConditions.atom_exponential_grid, 
            OneDSchrodingerEquationSolver.solver_exponential_grid ,
            l=i_orbi.l);
            #Update eigenvalue and eigenfunction in the basis set data structure.
            i_orbi.E=ei_temp;
            i_orbi.u=u_temp;
            E_eigen+= i_orbi.occu*ei_temp;
            #println(i_orbi.E)
            #println("--------------------------------")

        end
        #Update E_total_before from the E_total from the previous step.
        E_total_before= float(E_total);
        #Calculate density with new basis set.
        density_out= Density.calculate_density(basis);
        #Smooth the density with linear mixing (combination) of the previous and current densities.
        density_in= Density.linear_mixing(density_in, density_out, alpha=0.10);
        #Solve Poisson equation to find the new Hartree potential.
        V_hartree= OneDPoissonEquationSolver.solver_exponential_grid(Z, density_in, grid_stru);
        #Calculate new exchange and correlation potentials.
        V_x, E_xp, V_c, E_cp= ExchangeCorrelation.potentials(density_in);
        #Add exchange and correlation potentials.
        V_xcp= V_x .+ V_c;
        #Integrals to calculate energy components.

        V_xc= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_xcp)

        E_hartree= 0.5*4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_hartree);

        E_x= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, E_xp);

        E_c= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, E_cp);

        #Calculate total energy.
        E_total= E_eigen - E_hartree + E_x + E_c - V_xc;
        #println(E_total)
        #println("*****************************************")
    end
    @test abs(E_total - -37.425749) < 10.0e-5;
    # -37.425749 NIST value for LDA potential
end