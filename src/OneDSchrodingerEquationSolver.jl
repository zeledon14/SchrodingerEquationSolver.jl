using AutomaticDocstrings

module OneDSchrodingerEquationSolver

    using ..IntegralNumericalMethods
    using ..MathUtils



    function solver_uniform_integer_grid(E::Float64,u1::Float64,
        du1::Float64, i_1::Int64, u_end::Float64,
        du_end::Float64, i_end::Int64,
        v_effe::Vector{Float64},
        grid_stru::Any)::Tuple{Vector{Float64},Float64, Float64}

        grid_i::Vector{Float64}=grid_stru.grid_i;
        grid::Vector{Float64}=grid_stru.grid;
        dx_di::Vector{Float64}=grid_stru.dx_di;
        f::Vector{Float64}= 2.0.*(v_effe .- E);
        g=zeros(Float64, grid_stru.N);

        #find turn_pnts of of f, basically the clasical turning points of the effective density_potential
        #with restepect to the E proposed eigenvalue
        turn_pnts= MathUtils.indices_of_zeros_finder(f);
        if length(turn_pnts) ==  0
            throw(DomainError("the effective potential has no turning points 
            for the proposed energy eigenvalue, this means v_effe - E has no zeroes"));
        end
        #do forward integration of radial shcrodinger equation u
        u_fwd_temp= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_integer_grid(grid_i[i_1:end],
                g[i_1:end],f[i_1:end],
                dx_di[i_1:end],u1,du1);
        u_fwd= zeros(Float64, grid_stru.N);
        u_fwd[i_1:end].= u_fwd_temp;
        u_fwd_temp=Nothing;

        #do backward integreation of the radial shcrodinger equation u 
        u_bwd_temp= reverse(IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5_on_integer_grid(reverse(grid_i[1:i_end]),
        g,reverse(f[1:i_end]),
        reverse(dx_di[1:i_end]),u_end,du_end));
        u_bwd= zeros(Float64, grid_stru.N);
        u_bwd[1:i_end]= u_bwd_temp;
        u_bwd_temp=Nothing;

        #rescale u_fwd, u_bwd to make u_fwd[turn_pnts[1]] = u_bwd[turn_pnts[1]] = 1
        u_fwd, u_bwd= MathUtils.rescale_to_unity_at_turning_point!(u_fwd, u_bwd, turn_pnts[1]);
        #merge solutions
        u_merged, merge_value, merge_ratio= MathUtils.merge_solutions(u_fwd, u_bwd, grid, turn_pnts[1]);
        u_merged= MathUtils.normalize!(u_merged, grid);
        return u_merged, merge_value, merge_ratio
    end
end