using AutomaticDocstrings

module OneDSchrodingerEquationSolver

using ..IntegralNumericalMethods
using ..MathUtils

raw"""
    solver(E::Float64,init_valu1_fwrd::Float64,
    init_valu2_fwrd::Float64, init_valu1_bwrd::Float64,
    init_valu2_bwrd::Float64, v_effe::Vector{Float64},
    grid::Vector{Float64})::Tuple{Vector{Float64},Float64}

    Returns the normalized wave function (u(r)) of one dimensional shcrodinger equation of the form
    \ frac{-1}{2} \ frac{d^2 u(r)}{dr^2} + V_{effe} u(r) = E u(r) 
    For numerical stability and accuracy is integrated forward and backward and merged at the first 
    classical turning point. The return solution is normalized.
    It also returns the merging value to help in the finding of the energy eigenvalue.

**Inputs:**
- `E::Float64`: Value to use as energy eigenvalue
- `init_valu1_fwrd::Float64`: Value of u(grid[1]) to use as first initial value in forward solution.
- `init_valu2_fwrd::Float64`: Value of u(grid[2]) to use as second initial value in forward solution.
- `init_valu1_bwrd::Float64`: Value of u(grid[end]) to use as first initial value in backward solution.
- `init_valu2_bwrd::Float64`: Value of u(grid[end -1]) to use as second initial value in backward solution.
- `v_effe::Vector{Float64}`: Vector with the values of the effective potential over the grid.
- `grid::Vector{Float64}`: Vector with the grid values.
"""

function solver(E::Float64,init_valu1_fwrd::Float64,
    init_valu2_fwrd::Float64, init_valu1_bwrd::Float64,
    init_valu2_bwrd::Float64, v_effe::Vector{Float64},
    grid::Vector{Float64})::Tuple{Vector{Float64},Float64}

    f::Vector{Float64}= 2.0.*(v_effe .- E);
    g=zeros(Float64, size(f)[1]);

    #find turn_pnts of of f, basically the clasical turning points of the effective density_potential
    #with restepect to the E proposed eigenvalue
    turn_pnts= MathUtils.indices_of_zeros_finder(f);
    if length(turn_pnts) ==  0
        throw(DomainError("the effective potential has no turning points 
        for the proposed energy eigenvalue, this means v_effe - E has no zeroes"));
    end

    #do forward integration of radial shcrodinger equation u
    u_fwd= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5(grid,g,f,
    init_valu1_fwrd,init_valu2_fwrd);
    #do backward integreation of the radial shcrodinger equation u 
    u_bwd= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5(reverse(grid),g,reverse(f),
    init_valu1_bwrd,init_valu2_bwrd));
    #rescale u_fwd, u_bwd to make u_fwd[turn_pnts[1]] = u_bwd[turn_pnts[1]]
    u_fwd, u_bwd= MathUtils.rescale!(u_fwd, u_bwd, turn_pnts[1]);
    #merge solutions
    u_merged, merge_value= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
    u_merged= MathUtils.normalize!(u_merged, grid);
    return u_merged, merge_value
end

end