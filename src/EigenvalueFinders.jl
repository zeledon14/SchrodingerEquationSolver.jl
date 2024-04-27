
module EigenvalueFinders
using ..IntegralNumericalMethods
using ..MathUtils
import ..OneDSchrodingerEquationSolver as odses


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
function find_eigenvalue_intervals(energy_grid::Vector{Float64},v_effe::Vector{Float64}, grid::Vector{Float64}, 
    initial_condition_function::Function;
    l::Int64=0, integrador_type::String="RK4_PCABM5")::Tuple{Vector{Tuple{Float64,Float64}},Vector{Tuple{Float64,Float64}}}

    E_N= size(energy_grid)[1]
    merg_valu_of_E=zeros(Float64, E_N);
    for (i, ei) in enumerate(energy_grid)

        init_valu1_fwrd, init_valu2_fwrd,
        init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, ei, l);

        u_merged, merge_value= odses.solver(ei,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid, integrador_type);
        merg_valu_of_E[i]=merge_value;
    end
    ener_indx= MathUtils.indices_of_zeros_finder(merg_valu_of_E);
    #check that merg_valu_of_E is_continuous_enough around 


    good_intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx if MathUtils.is_continuous_enough(merg_valu_of_E, i, E_N)]
    bad_intervals=[(energy_grid[i], energy_grid[i+1]) for i in ener_indx if MathUtils.is_continuous_enough(merg_valu_of_E, i, E_N) == false]


    return good_intervals, bad_intervals


end

function illinois_eigenvalue_finder(E_interval::Tuple{Float64, Float64},
    v_effe::Vector{Float64}, grid::Vector{Float64}, 
    initial_condition_function::Function;
    l::Int64=0, 
    N_max::Int64=1000, tolerance::Float64=10.0e-12,
    integrador_type::String="RK4_PCABM5")::Tuple{Vector{Float64}, Float64}
    i=0
    Ec_befo=10.0e2
    Ea=E_interval[1]
    Eb=E_interval[2]
    Ec=0.0
    init_valu1_fwrd, init_valu2_fwrd,
    init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ea, l);
    _, u0a= odses.solver(Ea,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid,integrador_type);
    init_valu1_fwrd, init_valu2_fwrd,
    init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Eb, l);
    _, u0b= odses.solver(Eb,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid,integrador_type);
    while i < N_max
        Ec=(Ea*u0b -Eb*u0a)/(u0b - u0a)
        if abs(Ec-Ec_befo) < tolerance
            break
        end
        init_valu1_fwrd, init_valu2_fwrd,
        init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ec, l);
        _, u0c= odses.solver(Ec,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
            init_valu2_bwrd, v_effe, grid,integrador_type);
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
    init_valu1_fwrd, init_valu2_fwrd,
    init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, Ec, l);
    u, _= odses.solver(Ec,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid,integrador_type);
    return u, Ec
end

end