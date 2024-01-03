
module EigenvalueFinders
using ..IntegralNumericalMethods
using ..MathUtils
import ..OneDSchrodingerEquationSolver as odses



function find_eigenvalue_intervals(energy_grid::Vector{Float64},v_effe::Vector{Float64}, grid::Vector{Float64}, 
    initial_condition_function::Function, l::Int64=0)::Vector{Tuple{Float64,Float64}}

    E_N= size(energy_grid)[1]
    merg_valu_of_E=zeros(Float64, E_N);
    for (i, ei) in enumerate(energy_grid)

        init_valu1_fwrd, init_valu2_fwrd,
        init_valu1_bwrd, init_valu2_bwrd =initial_condition_function(grid, ei, l);

        u_merged, merge_value= odses.solver(ei,init_valu1_fwrd,init_valu2_fwrd, init_valu1_bwrd,
        init_valu2_bwrd, v_effe, grid);
        merg_valu_of_E[i]=merge_value;
    end
    ener_indx= MathUtils.indices_of_zeros_finder(merg_valu_of_E);
    out=[(energy_grid[i], energy_grid[i+1]) for (i) in ener_indx]
    return out


end

function illinois_eigenvalue_finder_from_guess(E_guess::Float64,
    E_plus::Float64, E_minu::Float64,
    v_effe::Vector{Float64}, grid::Vector{Float64}, 
    initial_condition_function::Function,
    N_max::Int64=300, tolerance::Float64=10.0e-10)
    i=0
    Ec_befo=E_guess
    Ea=nodes[1]
    Eb=nodes[2]
    Ec=0.0
    _, _, _, u0a, _= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Ea)
    _, _, _, u0b, _= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Eb)
    while i < N_max
        Ec=(Ea*u0b -Eb*u0a)/(u0b - u0a)
        if abs(Ec-Ec_befo) < tolerance
            break
        end
        _, _, _, u0c, _= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Ec)
        if Integer(sign(u0c)) == Integer(sign(u0a))
            Ea=float(Ec)
            u0a=float(u0c)
            u0b=0.5*u0b
        else
            Eb=float(Ec)
            u0b=float(u0c)
            u0a=0.5*u0a
        end
        Ec_befo=Ec
        i+=1
    end
    u, ub, uf, _, node= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Ec)
    u= Utils.normalize(grid, u)
    return u, ub, uf, node, Ec
end

end