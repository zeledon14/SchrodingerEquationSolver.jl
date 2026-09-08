using AutomaticDocstrings
using Test

using SchrodingerEquationSolver
using SchrodingerEquationSolver
using SchrodingerEquationSolver:   Potentials, MathUtils, Hydrogen, InitialConditions,
                                   IntegralNumericalMethods, OneDSchrodingerEquationSolver
@testset "H_TEST" begin
    n=1;
    Z=1;
    l=0;
    r_max= 25.0;
    E=-Z^2/(2*n^2);

    #grid_stru= Grids.UniformGrid(10.0e-3, 15.0, 5000);
    grid_stru= Grids.ExponentialGrid(r_max, Z);
    grid::Vector{Float64}=grid_stru.grid;
    grid_i::Vector{Float64}=grid_stru.grid_i;
    dx_di::Vector{Float64}=grid_stru.dx_di;
    grid_sqrt::Vector{Float64}=grid_stru.grid_sqrt;
    h_u_s1= Hydrogen.u_s1_hydrogen(grid_stru.grid);

    r_min=grid_stru.grid[1];
    u1= h_u_s1[1];
    w1= (h_u_s1[2] - h_u_s1[1])/(grid_stru.grid[2] - grid_stru.grid[1]);

    u_end= h_u_s1[end];
    w_end= (h_u_s1[end-1] - h_u_s1[end])/(grid_stru.grid[end-1] - grid_stru.grid[end]);

    #Z=1;
    v_effe =(
        Potentials.angular_potential(l, grid) .+
        Potentials.coulomb_potential(Z, grid));

    f::Vector{Float64}= 2.0.*(v_effe .- E);
    turn_pnts= MathUtils.indices_of_zeros_finder(f);
    g=zeros(Float64, size(f)[1]);
    u_fwd= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_integer_grid(grid_i,g,f,dx_di,
    u1,w1);

    u_bwd= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_integer_grid(reverse(grid_i),g,
    reverse(f),reverse(dx_di),
    u_end,w_end));
    u_fwd, u_bwd= MathUtils.rescale_abs!(u_fwd, u_bwd, turn_pnts[1]);
        #merge solutions
    u_merged, merge_value, merge_ratio= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
    error= ((h_u_s1 .- u_merged).^2.0).^0.5;
    max_error= maximum(error);
    average_error= sum(error)/length(error);

    @test max_error < 4e-7
    @test average_error < 2.5e-8

end

@testset "H_test_initial_conditions_solver_uniform_integer_grid" begin


    Z=1;
    r_max= 25.0;
    grid_stru= Grids.ExponentialGrid(r_max, Z);
    grid::Vector{Float64}=grid_stru.grid;
    h_u_s1= Hydrogen.u_s1_hydrogen(grid_stru.grid);
    h_u_s1_norm= MathUtils.normalize!(h_u_s1, grid);
    l=0;
    n=1;
    E=-Z^2/(2*n^2);
    v_effe =(
        Potentials.angular_potential(l, grid) .+
        Potentials.coulomb_potential(Z, grid));

    f::Vector{Float64}= 2.0.*(v_effe .- E);
    g=zeros(Float64, size(f)[1]);
    r_min=grid_stru.grid[1];
    u1, w1 = InitialConditions.atom_like_poly_at_r_min(r_min,l=l,E=E);
    u_end, w_end=InitialConditions.exponential_decay_at_r_ref(grid, lenght(grid), l=l, E=E);
    u_merged, merge_value, merge_ratio=OneDSchrodingerEquationSolver.solver_uniform_integer_grid(E,u1,w1, 
    u_end,w_end,v_effe,grid_stru)
    err= ((h_u_s1_norm .- u_merged).^2.0).^0.5;
    max_err= maximum(err);
    println("max error: ", max_err);
    average_error= sum(err)/length(err);
    println("average error: ", average_error);
    @test max_err < 1.4e-11
    @test average_error < 8.1e-14
    
end