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

function solver_uniform_grid(E::Float64,u1::Float64,
    du1::Float64, u_end::Float64,
    du_end::Float64, end_i::Int64,
    v_effe::Vector{Float64},
    grid_stru::Any)::Tuple{Vector{Float64},Float64, Float64}

    grid::Vector{Float64}=grid_stru.grid;
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
    u_fwd= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_direct_initial(grid,g,f,
    u1,du1);
    #do backward integreation of the radial shcrodinger equation u 
    u_bwd= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_direct_initial(reverse(grid),g,reverse(f),
    u_end,du_end));
    #rescale u_fwd, u_bwd to make u_fwd[turn_pnts[1]] = u_bwd[turn_pnts[1]]
    u_fwd, u_bwd= MathUtils.rescale!(u_fwd, u_bwd, turn_pnts[1]);
    #merge solutions
    u_merged, merge_value, merge_ratio= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
    u_merged= MathUtils.normalize!(u_merged, grid);
    return u_merged, merge_value, merge_ratio
end

#function solver_exponential_grid solver_v_return_u
function solver_exponential_grid(E::Float64,v1::Float64,
    dv1::Float64, v_end::Float64,
    dv_end::Float64, end_i::Int64,
    v_effe::Vector{Float64},
    grid_stru::Any)::Tuple{Vector{Float64},Float64, Float64}

    a::Float64=grid_stru.a;
    b::Float64=grid_stru.b;
    grid_i::Vector{Float64}=grid_stru.grid_i;
    f::Vector{Float64}= 2.0.*(v_effe .- E);
    fv::Vector{Float64}= ((a*b.*exp.(b.*grid_i)).^2).*f .+ 0.25*b^2;
    total_size=size(fv)[1];
    gv=zeros(Float64, total_size);

    #find turn_pnts of of f, basically the clasical turning points of the effective density_potential
    #with restepect to the E proposed eigenvalue
    turn_pnts= MathUtils.indices_of_zeros_finder(f);
    if length(turn_pnts) ==  0
        throw(DomainError("the effective potential has no turning points 
        for the proposed energy eigenvalue, this means v_effe - E has no zeroes"));
    end
  
    #do forward integration of radial shcrodinger equation u
    v_fwd= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_direct_initial(grid_i,gv,fv,
    v1,dv1);
    v_bwd= zeros(Float64, total_size);
    #do backward integreation of the radial shcrodinger equation u 
    v_bwd_short= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_direct_initial(reverse(grid_i[1:end_i]),gv[1:end_i],reverse(fv[1:end_i]),
    v_end, dv_end));
    v_bwd[1:end_i]= v_bwd_short;
    #rescale u_fwd, u_bwd to make u_fwd[turn_pnts[1]] = u_bwd[turn_pnts[1]]
    v_fwd, v_bwd= MathUtils.rescale!(v_fwd, v_bwd, turn_pnts[1]);
    #merge solutions
    v_merged, merge_value, merge_ratio= MathUtils.merge_solutions(v_fwd, v_bwd, grid_i, turn_pnts[1]);
    v_merged= MathUtils.normalize_v!(v_merged, grid_stru);
    u_merged= v_merged.*exp.((0.5*b).*grid_i);
    return u_merged, merge_value, merge_ratio
end
end