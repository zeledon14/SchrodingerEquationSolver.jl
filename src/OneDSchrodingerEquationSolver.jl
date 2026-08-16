using AutomaticDocstrings

module OneDSchrodingerEquationSolver

    using ..IntegralNumericalMethods
    using ..MathUtils



    function solver_uniform_integer_grid(E::Float64,u1::Float64,
        du1::Float64, u_end::Float64,
        du_end::Float64,
        v_effe::Vector{Float64},
        grid_stru::Any)::Tuple{Vector{Float64},Float64, Float64}

        grid_i::Vector{Float64}=grid_stru.grid_i;
        grid::Vector{Float64}=grid_stru.grid;
        dx_di::Vector{Float64}=grid_stru.dx_di;
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
        u_fwd= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_integer_grid(grid_i,g,f,
        dx_di,u1,du1);
        #do backward integreation of the radial shcrodinger equation u 
        u_bwd= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_integer_grid(reverse(grid_i),
        g,reverse(f),
        reverse(dx_di),u_end,du_end));
        #rescale u_fwd, u_bwd to make u_fwd[turn_pnts[1]] = u_bwd[turn_pnts[1]] = 1
        u_fwd, u_bwd= MathUtils.rescale_to_unity_at_turning_point!(u_fwd, u_bwd, turn_pnts[1]);
        #merge solutions
        u_merged, merge_value, merge_ratio= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
        u_merged= MathUtils.normalize!(u_merged, grid);
        return u_merged, merge_value, merge_ratio
    end
end