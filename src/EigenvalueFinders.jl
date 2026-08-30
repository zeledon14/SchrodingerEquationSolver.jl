
module EigenvalueFinders
using ..IntegralNumericalMethods
using ..MathUtils
import ..OneDSchrodingerEquationSolver.solver_uniform_integer_grid as solver


function find_eigenvalue_intervals(energy_grid::Vector{Float64},v_effe::Vector{Float64}, grid_stru::Any, 
    initial_condition_function_x_min::Function,
    initial_condition_function_x_max::Function,
    l::Int64=0)::Vector{Tuple{Float64,Float64}}#Tuple{Vector{Tuple{Float64,Float64}}, Vector{Float64}}#Tuple{Vector{Tuple{Float64,Float64}},Vector{Tuple{Float64,Float64}}}

    r_min=grid_stru.grid[1];
    r_max=grid_stru.grid[end];
    merge_value_list= zeros(length(energy_grid));
    for (i, E) in enumerate(energy_grid)
        u1, w1 = initial_condition_function_x_min(r_min,l=l, E=E);
        u_end, w_end=initial_condition_function_x_max(r_max, l=l,E=E);
        u_merged, merge_value, merge_ratio=solver(E, u1, w1, u_end, w_end, 
            v_effe, grid_stru);
        merge_value_list[i] = merge_value;
    end
    ener_indx= MathUtils.indices_of_zeros_finder(merge_value_list);
    #the energy indx should not be 1 
#clean the potential eigenvalue segments
    ener_indx_indicator= zeros(length(ener_indx));
    for (i,indx) in enumerate(ener_indx)
        delta_E= energy_grid[indx] - energy_grid[indx-1];
        delta_x= merge_value_list[indx] - merge_value_list[indx-1];
        log_slop= log10(abs(delta_x)/abs(delta_E));
        if log_slop < 1.1
            ener_indx_indicator[i]=1;
            #println("E= ", E_grid_stru.grid[indx],  " E_1= ", E_grid_stru.grid[indx-1], " merge_value= ", merge_value_list[indx]);
        end
    end    
    out_intervals::Vector{Tuple{Float64,Float64}}=[(0.0,0.0) for _ in 1:sum(ener_indx_indicator)];
    count=1;
    for (i,indx) in enumerate(ener_indx)
        if ener_indx_indicator[i]==1
            out_intervals[count]=(energy_grid[indx-1], energy_grid[indx]);
            count+=1;
            #println("E= ", E_grid_stru.grid[indx],  " E_1= ", E_grid_stru.grid[indx-1], " merge_value= ", merge_value_list[indx]);
        end
    end    
    return out_intervals#intervals, merge_ratio_of_E
end

function illinois_eigenvalue_finder(E_interval::Tuple{Float64, Float64},
    v_effe::Vector{Float64}, grid_stru::Any, 
    initial_condition_function_x_min::Function,
    initial_condition_function_x_max::Function,
    l::Int64=0, 
    N_max::Int64=1000, tolerance::Float64=10.0e-14)::Tuple{Vector{Float64}, Float64}
    i=0
    r_min=grid_stru.grid[1];
    r_max=grid_stru.grid[end];
    Ec_befo=10.0e2
    Ea=E_interval[1]
    Eb=E_interval[2]
    Ec=0.0
    u1, w1 = initial_condition_function_x_min(r_min,l=l, E=Ea);
    u_end, w_end=initial_condition_function_x_max(r_max, l=l,E=Ea);
    _, u0a, _=solver(Ea, u1, w1, u_end, w_end, 
        v_effe, grid_stru);
    #y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Ea, l);
    #_, u0a, _= solver(Ea, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
    #init_valu1_fwrd, init_valu2_fwrd,
    #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ea, l);
    #_, u0a= OneDSchrodingerEquationSolver.solver(Ea,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
    #    init_valu2_bwrd, v_effe, grid,integrador_type);
    u1, w1 = initial_condition_function_x_min(r_min,l=l, E=Eb);
    u_end, w_end=initial_condition_function_x_max(r_max, l=l,E=Eb);
    _, u0b, _=solver(Eb, u1, w1, u_end, w_end, 
        v_effe, grid_stru);
    #y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Eb, l);
    #_, u0b, _= solver(Eb, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
    #init_valu1_fwrd, init_valu2_fwrd,
    #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Eb, l);
    #_, u0b= OneDSchrodingerEquationSolver.solver(Eb,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
    #    init_valu2_bwrd, v_effe, grid,integrador_type);
    while i < N_max
        Ec=(Ea*u0b -Eb*u0a)/(u0b - u0a)
        if abs(Ec-Ec_befo) < tolerance
            break
        end
        u1, w1 = initial_condition_function_x_min(r_min,l=l, E=Ec);
        u_end, w_end=initial_condition_function_x_max(r_max, l=l,E=Ec);
        _, u0c, _=solver(Ec, u1, w1, u_end, w_end, 
            v_effe, grid_stru);        
        #y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Ec, l);
        #_, u0c,_= solver(Ec, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
        #init_valu1_fwrd, init_valu2_fwrd,
        #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ec, l);
        #_, u0c= OneDSchrodingerEquationSolver.solver(Ec,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        #    init_valu2_bwrd, v_effe, grid,integrador_type);
        if Integer(sign(u0c)) == Integer(sign(u0a))
            Ea=float(Ec)
            u0a=float(u0c)
            u0b=0.5*u0b
            #println(" 1  ", Ec)
        else
            Eb=float(Ec)
            u0b=float(u0c)
            u0a=0.5*u0a
            #println(" 2  ", Ec)
        end
        Ec_befo=Ec
        i+=1
    end
    u1, w1 = initial_condition_function_x_min(r_min,l=l, E=Ec);
    u_end, w_end=initial_condition_function_x_max(r_max, l=l,E=Ec);
    u, _, merge_ratio=solver(Ec, u1, w1, u_end, w_end, 
        v_effe, grid_stru);    
    #y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Ec, l);
    #u, _, merge_ratio= solver(Ec, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
    #init_valu1_fwrd, init_valu2_fwrd,
    #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ec, l);
    #u, _= OneDSchrodingerEquationSolver.solver(Ec,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
    #    init_valu2_bwrd, v_effe, grid,integrador_type);
    return u, Ec
end


function guess_energy_interval(eigen_before::Float64, V_effe_max::Float64, 
    V_effe_min::Float64, left_scale::Float64=0.35,
    right_scale::Float64=0.15)::Tuple{Float64,Float64}
    #TO DO   CHECK THAT THE INTERVAL HAS A SOFT EIGENVALUE
    E_guess_max= eigen_before - left_scale*eigen_before;
    E_guess_min= eigen_before + right_scale*eigen_before;
    while E_guess_max > V_effe_max
        if E_guess_max> 0
            E_guess_max= E_guess_max - 0.1*E_guess_max;
        else
            E_guess_max= E_guess_max + 0.1*E_guess_max;
        end
    end

    while E_guess_min < V_effe_min
        if E_guess_min < 0.0
            E_guess_min = E_guess_min - 0.1*E_guess_min;
        else
            E_guess_min = E_guess_min + 0.1*E_guess_min;
        end
    end
    return (E_guess_min, E_guess_max)
    
end

end