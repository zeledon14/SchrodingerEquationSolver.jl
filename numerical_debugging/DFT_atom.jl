
# Add src/ to the LOAD_PATH so Julia can find SchrodingerEquationSolver
#push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))
using SchrodingerEquationSolver
using SchrodingerEquationSolver: Grids, Potentials, MathUtils, Hydrogen, InitialConditions,
                                 OneDSchrodingerEquationSolver, OneDPoissonEquationSolver,
                                 EigenvalueFinders, AtomBasisSet, Density, ExchangeCorrelation
using Plots
using CSV
using DataFrames
using PrettyTables

#Define parameters and produce an exponential grid.
r_max::Float64=15;#Max radius of space grid.
Z::Int64=2; #Atomic number, also used as the charge of coulomb potential.

df_dft=DataFrame(CSV.File(joinpath("/home/arturo_hernandez/Downloads/nist_dft_tables/dft_Z_$Z/$Z.csv")));
# Extract first row as column names
col_names = Symbol.(collect(df_dft[1, :]));  # Convert strings to Symbols
# Remove the first row and assign new column names
df_dft = df_dft[2:end, :];
rename!(df_dft, col_names);
energy_target = Dict(zip(replace.(df_dft.Energy, " =" => "") , parse.(Float64, df_dft.LDA)));

#grid definition
grid_stru= Grids.init_exponential_grid_structure(r_max, Z);
N=grid_stru.N;
print("grid size ", N)

#Initialization of potentials and energies
let
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

E_kinetic::Float64= 0.0;
#Initializing basis set data structure
basis= AtomBasisSet.init_atom_basis_set(Z, grid_stru.grid);

#Energy minimization loop 
while abs(E_total - E_total_before) > 10.0e-8
    E_eigen=0.0;
    #Loop over every orbital to solve independent particle Schrodinger equation.
    for i_orbi in basis.orbitals
        println(i_orbi.name)
        #angular potential for l orbital
        V_angu= Potentials.angular_potential(i_orbi.l, grid_stru.grid);
        #Assemble effective potential.
        V_effe= V_colu .+ V_angu .+ V_hartree .+ V_x .+ V_c;
        V_effe_max= maximum(V_effe)
        V_effe_min= minimum(V_effe)
        println(V_effe_max);
        println(V_effe_min);
        energy_interval= EigenvalueFinders.guess_energy_interval(i_orbi.E, V_effe_max, V_effe_min);
        println("energy_interval = $energy_interval");
        E_grid= Grids.uniform_grid(energy_interval[1], energy_interval[2], 300); #List with the energy grid points.

        E_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, V_effe, grid_stru,
                    InitialConditions.atom_exponential_grid,
                        OneDSchrodingerEquationSolver.solver_exponential_grid, numb_inter=1);
        println("E_intervals = $E_intervals");
                        
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
    V_no_kinetic= V_colu .+ V_hartree .+ V_x .+ V_c;
    #Calculate density with new basis set.
    density_out= Density.calculate_density(basis);
    #Smooth the density with linear mixing (combination) of the previous and current densities.
    density_in= Density.linear_mixing(density_in, density_out, alpha=0.110);

    #Solve Poisson equation to find the new Hartree potential.
    V_hartree= OneDPoissonEquationSolver.solver_exponential_grid(Z, density_in, grid_stru);
    #Calculate new exchange and correlation potentials.
    V_x, E_xp, V_c, E_cp= ExchangeCorrelation.potentials(density_in);
    #Add exchange and correlation potentials.
    V_xcp= V_x .+ V_c;
    #Integrals to calculate energy components.
    E_kinetic= E_eigen - 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_no_kinetic);
    #V_xc= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_xcp)

    E_hartree= 0.5*4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_hartree);

    E_x= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, E_xp);

    E_c= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, E_cp);
    E_colu= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_colu);
    #Calculate total energy.
    E_total= E_kinetic + E_hartree + E_x + E_c + E_colu;
    #println(E_total)
    #println("*****************************************")

end
println("Ekin ", E_kinetic)
println("Ecoul ", E_hartree)
println("Exc ", (E_x + E_c))
println("Etot ", E_total)
pred_energy_dict= Dict("Ekin"=> E_kinetic, "Exc"=>(E_x + E_c),
"Ecoul"=>E_hartree, "Etot"=>E_total)
for i_orbi in basis.orbitals
    merge!(pred_energy_dict, Dict(i_orbi.name => i_orbi.E))
    println(i_orbi.name)
    println(i_orbi.E)
end
end#let