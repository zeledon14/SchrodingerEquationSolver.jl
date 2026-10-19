using SchrodingerEquationSolver
using SchrodingerEquationSolver: Grids, Potentials, MathUtils, Hydrogen, InitialConditions,
                                 OneDSchrodingerEquationSolver, OneDPoissonEquationSolver,
                                 EigenvalueFinders, BasisSets, Density, ExchangeCorrelation,
                                 PulayDensity
import ..OneDSchrodingerEquationSolver.solver_uniform_integer_grid as solver
using Plots
using Libxc
using CSV

let 
    n=3
    Z=92
    l=1
    r_max= 50.0
    E_0=-Z^2/(2*n^2)
    E_0= E_0 +(0.15)*E_0

    #grid_stru= Grids.UniformGrid(10.0e-3, 15.0, 5000)
    grid_stru= Grids.ExponentialGrid(r_max, Z)
    grid::Vector{Float64}=grid_stru.grid
    grid_i::Vector{Float64}=grid_stru.grid_i
    dx_di::Vector{Float64}=grid_stru.dx_di
    grid_sqrt::Vector{Float64}=grid_stru.grid_sqrt
    #basis= BasisSets.AtomBasisSet(Z, grid)
    initial_condition_function_x_min = InitialConditions.atom_like_poly_at_r_min
    initial_condition_function_x_max = InitialConditions.exponential_decay_at_r_ref


    v_effe =(
        Potentials.angular_potential(l, grid) .+
        Potentials.coulomb_potential(Z, grid))


    mvi=-10.0
    mv_lst=-20.0
    nodesi=-5.0
    node_lst=-6.0
    E_lst=-10.0
    E_nxt=-10.0
    Ei=E_0#initial guess
    target_nodes=n-l-1
    mv_tol=10.0e-15
    steps_max=10000
    E_dict=Dict([("max",Dict([("E",-10.0),("mv",-10.0)])), 
    ("min",Dict([("E",-10.0),("mv",-10.0)]))])
    node_search=true
    gradient_decent=true
    for i in 1:steps_max
        abs(mvi) < mv_tol && nodesi == target_nodes && break
    #while abs(mvi) > mv_tol && i < steps_max
        u1, w1, i_1 = initial_condition_function_x_min(grid, 1,l=l, E=Ei)
        u_end, w_end, i_end =initial_condition_function_x_max(grid, grid_stru.N, l=l,E=Ei)
        u_merged, mvi, nodesi=solver(Ei, u1, w1, i_1, u_end, w_end, i_end,
                                    v_effe, grid_stru)
        if i == 1
            mv_lst= mvi
            node_lst=nodesi
        end
        if nodesi != target_nodes && node_search 
            if nodesi < target_nodes
                E_nxt= Ei +0.1*abs(Ei)
            else
                E_nxt= Ei -0.1*abs(Ei)
            end
        else
            node_search && (mv_lst= mvi)
            if sign(mv_lst) == sign(mvi) && gradient_decent
                node_search=false
                E_nxt= Ei +0.05*sign(mvi)*abs(Ei)
                E_dict["max"]["E"]=Ei
                E_dict["max"]["mv"]=mvi
                E_dict["min"]["E"]=Ei
                E_dict["min"]["mv"]=mvi
            else #illinois algorithm
                gradient_decent=false
                if sign(E_dict["max"]["mv"]) == sign(mvi) 
                    E_dict["max"]["E"]=Ei
                    E_dict["max"]["mv"]=mvi
                    E_dict["min"]["mv"]=0.5*E_dict["min"]["mv"]
                else
                    E_dict["min"]["E"]=Ei
                    E_dict["min"]["mv"]=mvi
                    E_dict["max"]["mv"]=0.5*E_dict["max"]["mv"]
                end


                E_nxt= (E_dict["min"]["E"]*E_dict["max"]["mv"] - E_dict["max"]["E"]*E_dict["min"]["mv"])/(E_dict["max"]["mv"] - E_dict["min"]["mv"])
            end
        end
        mv_lst=mvi
        node_lst=nodesi
        E_lst=Ei
        Ei=E_nxt
        i+=1
    end

    u1, w1, i_1 = initial_condition_function_x_min(grid, 1,l=l, E=Ei)
    u_end, w_end, i_end =initial_condition_function_x_max(grid, grid_stru.N, l=l,E=Ei)
    u_merged, mvi, nodesi=solver(Ei, u1, w1, i_1, u_end, w_end, i_end,
                                v_effe, grid_stru)
    println(Ei)
    plot(grid_i, u_merged, title="$n , $l orbital")

end
#=
    E_grid= Grids.UniformGrid(1.03*E_0, 0.97*E_0, 250).grid
nodes_list=zeros(length(E_grid))
merge_value_list=zeros(length(E_grid))

for (i,E_i ) in enumerate(E_grid)
    merge_value_list[i] = merge_value
    nodes_list[i]= nodes
end

p1 = plot(E_grid, merge_value_list, title="merge_value_list vs energy_grid")
p2 = plot(E_grid, nodes_list, title="nodes_list vs energy_grid")
display(plot(p1, p2, layout=(2,1), size=(800,700)))
=#