
module EigenvalueFinders
using ..IntegralNumericalMethods
using ..MathUtils
import ..OneDSchrodingerEquationSolver as OneDSchrodingerEquationSolver


"""
    find_eigenvalue_intervals(energy_grid::Vector{Float64}, v_effe::Vector{Float64}, 
    grid::Vector{Float64}, initial_condition_function::Function, l::Int64 = 0)

Finds the intervals (E_{min}, E_{max}) in which an eigenvalue exist
if they exist in the given energy grid.

**Inputs:**
- `energy_grid::Vector{Float64}`: The grid of energy values to use to solve 
                                the equation as proposed eigenvalues.
- `v_effe::Vector{Float64}`: Effective potential to solve equation.
- `grid::Vector{Float64}`: The space grid where the functins are deffined.
- `initial_condition_function::Function`: Function to calculate initial values for 
                                        Schrodinger equation.
- `l::Int64`: Angular quantum number, deffault value 0
"""
function find_eigenvalue_intervals_old(energy_grid::Vector{Float64},v_effe::Vector{Float64}, grid::Vector{Float64}, 
    initial_condition_function::Function;
    l::Int64=0, integrador_type::String="RK4_PCABM5")::Tuple{Vector{Tuple{Float64,Float64}},Vector{Tuple{Float64,Float64}}}

    E_N= size(energy_grid)[1]
    merg_valu_of_E=zeros(Float64, E_N);
    for (i, ei) in enumerate(energy_grid)

        init_valu1_fwrd, init_valu2_fwrd,
        init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, ei, l);

        u_merged, merge_value= OneDSchrodingerEquationSolver.solver(ei,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid, integrador_type);
        merg_valu_of_E[i]=merge_value;
    end
    ener_indx= MathUtils.indices_of_zeros_finder(merg_valu_of_E);
    #check that merg_valu_of_E is_continuous_enough around 


    good_intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx if MathUtils.is_continuous_enough(merg_valu_of_E, i, E_N)]
    bad_intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx if MathUtils.is_continuous_enough(merg_valu_of_E, i, E_N) == false]


    return good_intervals, bad_intervals


end


function find_all_eigenvalue_intervals(energy_grid::Vector{Float64},v_effe::Vector{Float64}, grid_stru::Any, 
    initial_condition_function::Function,
    solver::Function;
    l::Int64=0)::Vector{Tuple{Float64,Float64}}#Tuple{Vector{Tuple{Float64,Float64}}, Vector{Float64}}#Tuple{Vector{Tuple{Float64,Float64}},Vector{Tuple{Float64,Float64}}}

    E_N= size(energy_grid)[1]
    merg_valu_of_E=zeros(Float64, E_N);
    merge_ratio_of_E=zeros(Float64, E_N);
    for (i, ei) in enumerate(energy_grid)

        v1, dv1, v_end, dv_end, end_i=initial_condition_function(grid_stru, ei, l);

        u_merged, merge_value, merge_ratio= solver(ei, v1, dv1, v_end, dv_end, end_i,v_effe, grid_stru);
        merg_valu_of_E[i]=merge_value;
        merge_ratio_of_E[i]=merge_ratio;
    end
    ener_indx= MathUtils.indices_of_zeros_finder(merg_valu_of_E);
    #check that merg_valu_of_E is_continuous_enough around 


    #good_intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx if MathUtils.is_continuous_enough(merg_valu_of_E, i, E_N)]
    #bad_intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx if MathUtils.is_continuous_enough(merg_valu_of_E, i, E_N) == false]
    intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx]
    merge_ratio_of_E=zeros(Float64, size(intervals)[1]);
    for (i, i_interval) in enumerate(intervals)
        ei= 0.5*(i_interval[1] + i_interval[2])
        v1, dv1, v_end, dv_end, end_i=initial_condition_function(grid_stru, ei, l);

        u_merged, merge_value, merge_ratio= solver(ei, v1, dv1, v_end, dv_end, end_i,v_effe, grid_stru);
        merge_ratio_of_E[i]=merge_ratio;
    end

    out_intervals=[(i_interval) for (i, i_interval) in enumerate(intervals) if merge_ratio_of_E[i] < 1.25];
    return out_intervals#intervals, merge_ratio_of_E


end

function find_eigenvalue_intervals(energy_grid::Vector{Float64},v_effe::Vector{Float64}, grid_stru::Any, 
    initial_condition_function::Function,
    solver::Function;
    l::Int64=0,
    numb_inter::Int64=0)::Vector{Tuple{Float64,Float64}}#Tuple{Vector{Tuple{Float64,Float64}}, Vector{Float64}}#Tuple{Vector{Tuple{Float64,Float64}},Vector{Tuple{Float64,Float64}}}

    if numb_inter == 0
        return find_all_eigenvalue_intervals(energy_grid,v_effe, grid_stru, initial_condition_function,solver;l);
    else
        out_intervals::Vector{Tuple{Float64,Float64}}=[(0.0,0.0) for _ in 1:numb_inter];
        intervals_count=1;

        E_N= size(energy_grid)[1]

        e_befo=energy_grid[1];
        #println("e_befo  ", e_befo)
        v1, dv1, v_end, dv_end, end_i=initial_condition_function(grid_stru, e_befo , l);

        u_merged, merge_value_befo, merge_ratio_befo= solver(e_befo, v1, dv1, v_end, dv_end, end_i,v_effe, grid_stru);

        i=2;
        while i <=E_N && intervals_count <= numb_inter
            ei=energy_grid[i];
            #println("valiu ei ", ei)
            v1, dv1, v_end, dv_end, end_i=initial_condition_function(grid_stru, ei, l);

            u_merged, merge_value_curr, merge_ratio_curr= solver(ei, v1, dv1, v_end, dv_end, end_i,v_effe, grid_stru);
            if Int(sign(merge_value_befo)) != Int(sign(merge_value_curr))
                ei= 0.5*(e_befo + ei)
                #println("ratio ei ", ei)
                v1, dv1, v_end, dv_end, end_i=initial_condition_function(grid_stru, ei, l);
                u_merged, merge_value, merge_ratio= solver(ei, v1, dv1, v_end, dv_end, end_i,v_effe, grid_stru);
                #if merge_ratio < 1.25
                if merge_ratio < 1.1
                    out_intervals[intervals_count]=(e_befo,ei);
                    intervals_count+=1
                    #print("intervals_count ", intervals_count)
                end
            end
            e_befo= float(energy_grid[i]);
            i+=1;
        end
        return out_intervals#intervals, merge_ratio_of_E
    end


end

function illinois_eigenvalue_finder(E_interval::Tuple{Float64, Float64},
    v_effe::Vector{Float64}, grid_stru::Any, 
    initial_condition_function::Function,
    solver::Function;
    l::Int64=0, 
    N_max::Int64=1000, tolerance::Float64=10.0e-12)::Tuple{Vector{Float64}, Float64}
    i=0
    Ec_befo=10.0e2
    Ea=E_interval[1]
    Eb=E_interval[2]
    Ec=0.0
    y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Ea, l);
    _, u0a, _= solver(Ea, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
    #init_valu1_fwrd, init_valu2_fwrd,
    #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ea, l);
    #_, u0a= OneDSchrodingerEquationSolver.solver(Ea,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
    #    init_valu2_bwrd, v_effe, grid,integrador_type);
    y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Eb, l);
    _, u0b, _= solver(Eb, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
    #init_valu1_fwrd, init_valu2_fwrd,
    #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Eb, l);
    #_, u0b= OneDSchrodingerEquationSolver.solver(Eb,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
    #    init_valu2_bwrd, v_effe, grid,integrador_type);
    while i < N_max
        Ec=(Ea*u0b -Eb*u0a)/(u0b - u0a)
        if abs(Ec-Ec_befo) < tolerance
            break
        end
        y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Ec, l);
        _, u0c,_= solver(Ec, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
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
    y0_0, y1_0, y0_end, y1_end, end_i=initial_condition_function(grid_stru, Ec, l);
    u, _, merge_ratio= solver(Ec, y0_0, y1_0, y0_end, y1_end, end_i, v_effe, grid_stru);
    #init_valu1_fwrd, init_valu2_fwrd,
    #init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ec, l);
    #u, _= OneDSchrodingerEquationSolver.solver(Ec,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
    #    init_valu2_bwrd, v_effe, grid,integrador_type);
    return u, Ec
end


function guess_energy_interval(eigen_before::Float64, V_effe_max::Float64, 
    V_effe_min::Float64, left_scale::Float64=0.30,
    right_scale::Float64=0.01)::Tuple{Float64,Float64}
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