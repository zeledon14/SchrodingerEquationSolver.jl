
module EigenvalueFinders
using ..IntegralNumericalMethods
using ..MathUtils
import ..OneDSchrodingerEquationSolver as odses


function illinois_eigenvalue_finder_from_guess(E_guess::Float32,
    E_plus::Float32, E_minu::Float32,
    v_effe::Vector{Float32}, grid::Vector{Float32}, 
    initial_condition_function::Function,
    N_max::Int32=300, tolerance::Float32=10.0e-10)
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