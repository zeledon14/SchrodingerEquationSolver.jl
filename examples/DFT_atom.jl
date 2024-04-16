using Pkg;
Pkg.activate("../../SchrodingerEquationSolver");
include("../src/SchrodingerEquationSolver.jl");
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
import .ses.ExchangeCorrelation as ExchangeCorrelation

#Define parameters and produce an exponential grid.
r_max::Float64=50.0;#Max radius of space grid.
Z::Int64=6;#Atomic number, also used as the charge of coulomb potential.
b::Float64=0.001504; #By decreasing b the number of points in the grid increase.
grid= ses.Grids.exponential_grid(r_max, Z, b);
grid_sqrt= grid.^2.0; 
N= size(grid)[1]; #Number of points in the grid.

#Initialization of potentials and energies

#Initializing coulomb potential due to nuclei charge.
V_colu::Vector{Float64}= Potentials.coulomb_potential(Z, grid);
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
basis= AtomBasisSet.init_atom_basis_set(Z, grid);

#Energy minimization loop 
while abs(E_total - E_total_before) > 10.0e-12
    E_eigen=0.0;
    #Loop over every orbital to solve independent particle Schrodinger equation.
    for i_orbi in basis.orbitals
        #angular potential for l orbital
        V_angu= Potentials.angular_potential(i_orbi.l, grid);
        #Assemble effective potential.
        V_effe= V_colu .+ V_angu .+ V_hartree .+ V_x .+ V_c;
        #Energy grid to search for new eigenvalue, search around the previous eigenvalue.
        E_grid= Grids.uniform_grid(i_orbi.E - 0.5*i_orbi.E, i_orbi.E + 0.5*i_orbi.E, 150);
        #Find the energy intervals with potential eigenvalues.
        E_intervals, bad_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, V_effe,
                                                        grid,InitialConditions.atom, i_orbi.l);
        #Search for a solution to the independent particle Schrodinger equation for the first energy interval, 
        #there should be only one interval. Any other interval would be of a higher energy.
        u_temp, ei_temp= EigenvalueFinders.illinois_eigenvalue_finder(E_intervals[1], V_effe, grid,InitialConditions.atom);
        #Update eigenvalue and eigenfunction in the basis set data structure.
        i_orbi.E=ei_temp;
        i_orbi.u=u_temp;
        E_eigen+= i_orbi.occu*ei_temp;

    end
    #Update E_total_before from the E_total from the previous step.
    E_total_before= float(E_total);
    #Calculate density with new basis set.
    density_out= Density.calculate_density(basis);
    #Smooth the density with linear mixing (combination) of the previous and current densities.
    density_in= Density.linear_mixing(density_in, density_out, alpha=0.25);
    #Solve Poisson equation to find the new Hartree potential.
    V_hartree= odpes.solver(Z, density_in, grid);
    #Calculate new exchange and correlation potentials.
    V_x, E_xp, V_c, E_cp= ExchangeCorrelation.potentials(density_in);
    #Add exchange and correlation potentials.
    V_xcp= V_x .+ V_c;
    #Integrals to calculate energy components.
    V_xc= (4.0*pi)*MathUtils.integral((V_xcp.*density_in.*grid_sqrt), grid);
    E_hartree= 0.5*(4.0*pi)*(MathUtils.integral((V_hartree.*density_in.*grid_sqrt), grid));
    E_x= (4.0*pi)*MathUtils.integral((E_xp.*density_in.*grid_sqrt), grid);
    E_c= (4.0*pi)*MathUtils.integral((E_cp.*density_in.*grid_sqrt), grid);
    #C_in= (4.0*pi)*MathUtils.integral((density_in.*grid_sqrt), grid)
    #C_out= (4.0*pi)*MathUtils.integral((density_out.*grid_sqrt), grid)
    #Calculate total energy.
    E_total= E_eigen - E_hartree + E_x + E_c - V_xc;
end

println("total energy ", E_total)

save_path="../save_basis_set/test_write_basis.json";
AtomBasisSet.save_basis_set(basis, save_path);