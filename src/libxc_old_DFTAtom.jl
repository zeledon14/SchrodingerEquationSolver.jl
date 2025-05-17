
#using Pkg;
#Pkg.activate("../.");
module libxc_DFTAtom
using SchrodingerEquationSolver
using SchrodingerEquationSolver: Grids, Potentials, MathUtils, Hydrogen, InitialConditions,
                                 OneDSchrodingerEquationSolver, OneDPoissonEquationSolver,
                                 EigenvalueFinders, AtomBasisSet, Density, ExchangeCorrelation
using Plots
using CSV
using DataFrames
using PrettyTables
using PDFmerger
using LaTeXStrings
using Downloads  # for file deletion
using Printf
using Libxc

#module DFTAtom

    function calculate_atomic_basis_set(Z::Int64; r_max::Float64=50.0,
        potential_type::String="Free_atom", s::Float64= 200.0,
        r_onset::Float64= 4.00)::AtomBasisSet.atom_basis_set
        #Define parameters and produce an exponential grid.
        #r_max::Float64=50;#Max radius of space grid.
        #Z::Int64=8; #Atomic number, also used as the charge of coulomb potential.

        df_dft=DataFrame(CSV.File(joinpath(dirname(@__FILE__),"../nist_dft_tables/dft_Z_$Z/$Z.csv")));
        pdf_output=joinpath(dirname(@__FILE__),"../pdf_outputs/$(Z)_$(potential_type).pdf")
        # Extract first row as column names
        col_names = Symbol.(collect(df_dft[1, :]));  # Convert strings to Symbols
        # Remove the first row and assign new column names
        df_dft = df_dft[2:end, :];
        rename!(df_dft, col_names);
        energy_target = Dict(zip(replace.(df_dft.Energy, " =" => "") , parse.(Float64, df_dft.LDA)));

        #grid definition
        grid_stru= Grids.init_exponential_grid_structure(r_max, Z);
        N=grid_stru.N;
        xc = Functional(:lda_c_vwn);

        let
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

            E_kinetic::Float64= 0.0;

            V_conf::Vector{Float64}= zeros(Float64, N);
            if potential_type == "Confined_atom"
                #s= 200.0;#In hartree
                #r_onset= 4.00;#bohr radii
                V_conf= Potentials.Blum_confinement_potential(s, r_onset,grid_stru.grid);
            end
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
                    #V_effe= V_colu .+ V_angu .+ V_hartree .+ V_x .+ V_c.+ V_conf;
                    V_effe= V_colu .+ V_angu .+ V_hartree .+ V_conf .+ V_xcp;
                    V_effe_max= maximum(V_effe)
                    V_effe_min= minimum(V_effe)
                    #println(V_effe_max);
                    #println(V_effe_min);
                    energy_interval= EigenvalueFinders.guess_energy_interval(i_orbi.E, V_effe_max, V_effe_min);
                    #println(E_intervals);
                    E_grid= Grids.uniform_grid(energy_interval[1], energy_interval[2], 300); #List with the energy grid points.

                    E_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, V_effe, grid_stru,
                                InitialConditions.atom_exponential_grid,
                                    OneDSchrodingerEquationSolver.solver_exponential_grid, numb_inter=1);
                    #println(E_intervals);
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
                V_no_kinetic= V_colu .+ V_hartree .+ V_xcp;#V_x .+ V_c;
                #Calculate density with new basis set.
                density_out= Density.calculate_density(basis);
                #Smooth the density with linear mixing (combination) of the previous and current densities.
                density_in= Density.linear_mixing(density_in, density_out, alpha=0.110);

                #Solve Poisson equation to find the new Hartree potential.
                V_hartree= OneDPoissonEquationSolver.solver_exponential_grid(Z, density_in, grid_stru);
                #Calculate new exchange and correlation potentials.
                xc_all = evaluate(xc, rho=density_in);
                #V_x, E_xp, V_c, E_cp= ExchangeCorrelation.potentials(density_in);
                #Add exchange and correlation potentials.
                V_xcp= vec(xc_all.vrho);#V_x .+ V_c;
                E_xp= vec(xc_all.zk);
                #Integrals to calculate energy components.
                E_kinetic= E_eigen - 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_no_kinetic);
                #V_xc= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_xcp)

                E_hartree= 0.5*4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_hartree);

                E_x= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, E_xp);

                #E_c= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, E_cp);
                E_colu= 4.0*pi*MathUtils.energy_integral_exponential_grid(grid_stru, density_in, V_colu);
                #Calculate total energy.
                #E_total= E_kinetic + E_hartree + E_x + E_c + E_colu;
                E_total= E_kinetic + E_hartree + E_x + E_colu;

                #println(E_total)
                #println("*****************************************")

            end
            #println("Ekin ", E_kinetic)
            #println("Ecoul ", E_hartree)
            #println("Exc ", (E_x + E_c))
            #println("Etot ", E_total)
            pred_energy_dict= Dict("Ekin"=> E_kinetic, "Exc"=>(E_x + E_c),
            "Ecoul"=>E_hartree, "Etot"=>E_total)
            for i_orbi in basis.orbitals
                merge!(pred_energy_dict, Dict(i_orbi.name => i_orbi.E))
                #println(i_orbi.name)
                #println(i_orbi.E)
            end

            # Initialize arrays
            names = String[]
            E_pred = Float64[]
            E_targ = Float64[]

            for (key, value) in pred_energy_dict
                push!(names, key)
                push!(E_pred, value)
                push!(E_targ, energy_target[key])
            end

            # Compute absolute differences and percentage error
            E_diff = abs.(E_pred .- E_targ)
            Perc_error = E_diff ./ E_targ

            # Create table rows as vectors of strings
            atom_header = ["$(basis.Name)", "E in Hartree", "", "", ""]
            header = ["Name", "E NIST", "E pred", "|ΔE|", "|ΔE|/E NIST"]
            data_rows = [atom_header,header]
            for (n, targ, pred, diff, perc) in zip(names, E_targ, E_pred, E_diff, Perc_error)
                # Format each value nicely:
                push!(data_rows, [ n,
                                    @sprintf("%10.6f", targ),
                                    @sprintf("%10.6f", pred),
                                    @sprintf("%8.6f", diff),
                                    @sprintf("%12.6f", perc)
                                ])
            end

            # ----------------------------
            # Define table layout parameters
            # ----------------------------
            nrows = length(data_rows)     # number of rows (including header)
            ncols = length(header)        # number of columns

            # Define the horizontal boundaries for columns (you can adjust these widths)
            # For example, here we set 5 columns with different widths in arbitrary units.
            col_widths = [1.5, 2.5, 2.5, 2, 3]  # relative widths for each column
            # Compute cumulative positions
            col_x = [0.0]
            for w in col_widths
                push!(col_x, col_x[end] + w)
            end
            total_width = col_x[end]

            # Each row will have a fixed height:
            row_height = 1.0
            total_height = nrows * row_height

            # ----------------------------
            # Create the plot that will serve as the table "page"
            # ----------------------------
            p_table = plot(xlim=(0, total_width), ylim=(0, total_height),
                framestyle = :none,
                legend = false,
                grid = false,
                xticks = false,
                yticks = false,
                size = (800, Int(50 * nrows))  # adjust size as needed
            )
            # Add the title at the top center


            # Draw horizontal lines (rows)
            for i in 0:nrows
                plot!([0, total_width], [i * row_height, i * row_height], color = :black, lw = 1)
            end

            # Draw vertical lines (columns)
            for x in col_x
                plot!([x, x], [0, total_height], color = :black, lw = 1)
            end

            # Annotate each cell with the cell text (we assume center alignment in each cell)
            for (i, row) in enumerate(reverse(data_rows))  # reverse so that first row is at the top
                y_center = (i - 0.5) * row_height  # center of the row cell
                for j in 1:ncols
                    # Calculate the x-center for column j (between col_x[j] and col_x[j+1])
                    x_center = (col_x[j] + col_x[j+1]) / 2
                    text_str = row[j]
                    annotate!(p_table, x_center, y_center, text(text_str, 12, :center))
                end
            end

            # Save the table as a PDF
            table_page_file="table_page.pdf"
            savefig(p_table, table_page_file)

            orbital_files = String[]

            for i_orbi in basis.orbitals
                orbi_eigen = @sprintf("%.4f", i_orbi.E)
                plot_str = "Orbital $(i_orbi.name) eigen energy $orbi_eigen"

                p_orbi = plot(
                    grid_stru.grid_i, i_orbi.u,
                    title = plot_str,
                    xlabel = "i", ylabel = "u(r(i))",
                    legend = false
                )

                orbital_file = "orbital_$(i_orbi.name).pdf"
                savefig(p_orbi, orbital_file)
                push!(orbital_files, orbital_file)
            end



            all_files = [table_page_file; orbital_files...]
            final_pdf = pdf_output#"file.pdf"  # Final output file
            merge_pdfs(all_files, final_pdf)



            for f in all_files
                try
                    rm(f)
                catch e
                    @warn "Could not delete file $f" exception=(e, catch_backtrace())
                end
            end

            #println("✅ Successfully created $final_pdf with table and orbital plots!")
            basis_save_path=joinpath(dirname(@__FILE__),"../save_basis_set/$(potential_type)_z_$(Z)_r_max_$r_max.json")
            #joinpath(dirname(@__FILE__),"../save_basis_set/free_atom_z_$(Z)_r_max_$r_max.json")
            AtomBasisSet.save_basis_set(basis,basis_save_path)
        return basis
        end#let
    end

end

