using AutomaticDocstrings
using Test

using SchrodingerEquationSolver
using SchrodingerEquationSolver
using SchrodingerEquationSolver:   Potentials, MathUtils, Hydrogen, InitialConditions,
                                   IntegralNumericalMethods, OneDSchrodingerEquationSolver
@testset "CTest" begin
n=1;
Z=1;
l=0;
r_max= 25.0;
E=-Z^2/(2*n^2);

#grid_stru= Grids.init_uniform_grid_structure(10.0e-3, 15.0, 5000);
grid_stru= Grids.init_exponential_grid_structure(r_max, Z);
grid::Vector{Float64}=grid_stru.grid;
grid_i::Vector{Float64}=grid_stru.grid_i;
dr_di::Vector{Float64}=grid_stru.dr_di;
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
u_fwd_i= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_grid_i(grid_i,g,f,dr_di,
u1,w1);

u_bwd_i= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_grid_i(reverse(grid_i),g,
reverse(f),reverse(dr_di),
u_end,w_end));
u_fwd, u_bwd= MathUtils.rescale!(u_fwd, u_bwd, turn_pnts[1]);
    #merge solutions
u_merged, merge_value, merge_ratio= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
error= ((h_u_s1 .- u_merged).^2.0).^0.5;
max_error= maximum(error);
average_error= sum(error)/length(error);

@test max_error < 4e-7
@test average_error < 2.5e-8

end