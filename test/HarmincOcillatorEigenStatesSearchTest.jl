using AutomaticDocstrings
using Test

import .SchrodingerEquationSolver.Grids as Grids
import .SchrodingerEquationSolver.Potentials as Potentials
import .SchrodingerEquationSolver.Hydrogen as Hydrogen
import .SchrodingerEquationSolver.IntegralNumericalMethods as IntegralNumericalMethods
import .SchrodingerEquationSolver.MathUtils as MathUtils
import .SchrodingerEquationSolver.EigenvalueFinders as EigenvalueFinders

@testset "HarmincOcillatorEigenStatesSearchTest" begin
    #initial parameters
    r_min::Float64=-6.0;
    r_max::Float64=6.0;
    N=10000;
    w::Float64=1;
    m::Float64=1;
    #space grid creation
    grid= Grids.uniform_grid(r_min, r_max, N);
    #target values to compare calculations
    eigen_list_target::Vector{Float64}=[elem + 0.5 for elem in 0.0:3.0];
    u_wave_functions_target::Vector{Vector{Float64}}= [zeros(N) for _ in 1:4];
    harm_osci0::Vector{Float64}= exp.(-0.5.*abs.(grid).^2);
    u_wave_functions_target[1]= MathUtils.normalize!(harm_osci0, grid);

    harm_osci1::Vector{Float64}= 2*grid.*harm_osci0;
    u_wave_functions_target[2]= (-1.0).*MathUtils.normalize!(harm_osci1, grid);

    harm_osci2::Vector{Float64}= (4*grid.^2 .- 2).*harm_osci0;
    u_wave_functions_target[3]= MathUtils.normalize!(harm_osci2, grid);

    harm_osci3::Vector{Float64}= (8*grid.^3 .- 12*grid).*harm_osci0;
    u_wave_functions_target[4]= (-1.0).*MathUtils.normalize!(harm_osci3, grid);

    #define potential
    v_effe= Potentials.harmoic_oscilator_potential(w,m,grid);
    #define energy grid for eigenvalue search
    E_max::Float64=3.7;
    E_min::Float64=0.4;
    E_N::Int64=200;
    E_grid= Grids.uniform_grid(E_min, E_max, E_N);
    #serch interval of eigenvalues
    E_intervals, bad_intervals= EigenvalueFinders.find_eigenvalue_intervals(E_grid, v_effe, grid,InitialConditions.atom);
    numb_solu::Int64= size(E_intervals)[1];
    #test that the number of eigensates find is 4
    @test numb_solu == 4
    #define structures to calculated values of eigen energies
    #and eigen functions
    eigen_list::Vector{Float64}=zeros(numb_solu);
    u_wave_functions::Vector{Vector{Float64}}= [zeros(N) for _ in 1:numb_solu];
    for (i, ei_interval) in enumerate(E_intervals)
        u_temp, ei_temp= EigenvalueFinders.illinois_eigenvalue_finder(ei_interval, v_effe, grid,InitialConditions.atom);
        u_wave_functions[i]=u_temp;
        eigen_list[i]= ei_temp;
    end
    #test eigen energies
    for (i,ei) in enumerate(eigen_list_target)
        @test (abs(ei - eigen_list[i])) < 10.0e-9
    end
    #test eigen functions
    for (i, ui_wf) in enumerate(u_wave_functions_target)
        @test (MathUtils.error_difference(u_wave_functions[i],ui_wf)) < 10.0e-7
    end

end