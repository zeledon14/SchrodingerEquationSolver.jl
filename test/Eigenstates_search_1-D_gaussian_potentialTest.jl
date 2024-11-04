using AutomaticDocstrings
using Test

using SchrodingerEquationSolver
using SchrodingerEquationSolver: Grids, Potentials, MathUtils, Hydrogen, InitialConditions,
                                 OneDSchrodingerEquationSolver, OneDPoissonEquationSolver,
                                 EigenvalueFinders, AtomBasisSet, Density, ExchangeCorrelation
@testset "Eigenstates_search_1-D_gaussian_potentialTest" begin

# Space grid definition and creation
r_min::Float64=-12.0; #Where the space grid starts.
r_max::Float64=12.0; #Where the space grid ends.
N=40000; #Number of points in the space grind.
grid_stru= Grids.init_uniform_grid_structure(r_min, r_max, N); #Grid creation, grid is the list with the grid points.

v0::Float64=-10.0;

v_effe= Potentials.gaussian_potential(v0,grid_stru.grid); #the list with the gaussian potential

#Energy grid definition and creation. The system searches the energy grid for potential 
#values for the eigenvalues of the Schrodinger equation 
E_max::Float64=-0.05; #Maximal energy in the energy grid.
E_min::Float64=-9.98; #Minimal energy in the energy grid.
E_N::Int64=400; #Number of points in the energy grid.
E_grid= Grids.uniform_grid(E_min, E_max, E_N); #List with the energy grid points.

E_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, v_effe, grid_stru,
    InitialConditions.harmoic_oscillator_1D,
    OneDSchrodingerEquationSolver.solver_uniform_grid);

numb_solu::Int64= size(E_intervals)[1]; #Number of potential solutions in the energy grid.
eigen_list::Vector{Float64}=zeros(numb_solu); #Initializing the list that is going to hold the energy eigenvalues.
u_wave_functions::Vector{Vector{Float64}}= [zeros(N) for _ in 1:numb_solu]; #Initializing the list that holds the eigenfunctions list.

#Using Illinois algorithm to find the actual energy eigenvalue and eigenfunction for everyone of the energy intervals.
for (i, ei_interval) in enumerate(E_intervals)
    u_temp, ei_temp= EigenvalueFinders.illinois_eigenvalue_finder(ei_interval, v_effe, 
    grid_stru,InitialConditions.harmoic_oscillator_1D, 
    OneDSchrodingerEquationSolver.solver_uniform_grid);
    u_wave_functions[i]=u_temp;
    eigen_list[i]= ei_temp;
end

    @test abs(eigen_list[1] - -7.952670248742307) < 10.0e-6;
    @test abs(eigen_list[2] - -4.28060217561834) < 10.0e-6;
    @test abs(eigen_list[3] - -1.541950723244179) < 10.0e-6;
    @test abs(eigen_list[4] - -0.061136688426465795) < 10.0e-6;
    
end