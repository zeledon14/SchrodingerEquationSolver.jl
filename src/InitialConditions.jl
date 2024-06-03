using AutomaticDocstrings

module InitialConditions

"""
    atom(grid::Vector{Float64}, E::Float64, l::Int64)

Return the initial conditions for the electron u wave function for a hydrogenic atom
the initial conditions are based on asymptotic behavior for r-> 0 for the forward 
and r-> infinitive for backward.

**Inputs:**
- `grid::Vector{Float64}`: DESCRIPTION
- `E::Float64`: DESCRIPTION
- `l::Int64`: DESCRIPTION
"""
function atom(grid::Vector{Float64}, E::Float64=-0.5, 
              l::Int64=0)::Tuple{Float64,Float64,Float64,Float64}
    init_valu1_fwrd::Float64=grid[1]^(l+1.0)#u_s1_hydr_norm[1];
    init_valu2_fwrd::Float64=grid[2]^(l+1.0) #u_s1_hydr_norm[2];
    if sign(E) < 0.0
        lambda= (-2.0*E)^0.5;
    else
        lambda= (2.0*E)^0.5;
    end
    init_valu1_bwrd::Float64=grid[end]*exp(-1.0*lambda*grid[end]);#u_s1_hydr_norm[end];
    init_valu2_bwrd::Float64=grid[end-1]*exp(-1.0*lambda*grid[end-1]);#u_s1_hydr_norm[end-1];          
    return init_valu1_fwrd, init_valu2_fwrd, init_valu1_bwrd, init_valu2_bwrd
end

function u_r_near_0(grid::Vector{Float64}, i::Int64,l::Int64)::Float64
    return grid[i]^(l+1);
end

function du_dr_near_0(grid::Vector{Float64}, i::Int64,l::Int64)::Float64
    return (l+1)*grid[i]^l
end

function v_i_near_0(grid::Vector{Float64}, i::Int64,
    l::Int64, b::Float64)::Float64
    return u_r_near_0(grid, i, l)*exp(-0.5*b*i)
end

function dv_di_near_0(grid::Vector{Float64}, i::Int64,
    l::Int64, a::Float64,b::Float64)::Float64
    temp= du_dr_near_0(grid,i,l)*a*b*exp(0.5*b*i);
    temp1= 0.5*b*u_r_near_0(grid,i,l)*exp(-0.5*b*i);
    return temp - temp1
end

function u_r_end(grid::Vector{Float64}, E::Float64)::Float64
    if sign(E) < 0.0
        lambda= (-2.0*E)^0.5;
    else
        lambda= (2.0*E)^0.5;
    end

    return exp(-1.0*lambda*grid[end]);
end

function du_dr_end(grid::Vector{Float64}, E::Float64)::Float64
    if sign(E) < 0.0
        lambda= (-2.0*E)^0.5;
    else
        lambda= (2.0*E)^0.5;
    end

    return -1.0*lambda*exp(-1.0*lambda*grid[end]);
end


function v_i_end(grid::Vector{Float64}, E::Float64, b::Float64, i::Int64)::Float64
    return u_r_end(grid, E)*exp(-0.5*b*i)
end

function dv_di_end(grid::Vector{Float64},
    E::Float64, a::Float64,b::Float64, i::Int64)::Float64
    temp= du_dr_end(grid,E)*a*b*exp(0.5*b*i);
    temp1= 0.5*b*u_r_end(grid, E)*exp(-0.5*b*i);
    return temp - temp1
end

function atom_v(grid_stru::Any, E::Float64=-0.5, 
    l::Int64=0)::Tuple{Float64,Float64,Float64,Float64}
    end_i= size(grid_stru.grid_i)[1];
    v1= InitialConditions.v_i_near_0(grid_stru.grid, 1, l, grid_stru.b);
    dv1= InitialConditions.dv_di_near_0(grid_stru.grid,1, l, grid_stru.a, grid_stru.b);
    v_end= InitialConditions.v_i_end(grid_stru.grid, E, grid_stru.b,end_i);
    dv_end= InitialConditions.dv_di_end(grid_stru.grid, E, grid_stru.a, grid_stru.b, end_i);         
return v1, dv1, v_end, dv_end
end

function harmoic_oscillator(grid::Vector{Float64},  E::Float64=-0.5, 
    l::Int64=0)::Tuple{Float64,Float64,Float64,Float64}
    init_valu1_fwrd::Float64=exp(-0.5*abs(grid[1])^2);
    init_valu2_fwrd::Float64=exp(-0.5*abs(grid[2])^2);
    init_valu1_bwrd::Float64=exp(-0.5*abs(grid[end])^2);
    init_valu2_bwrd::Float64=exp(-0.5*abs(grid[end-1])^2);
    return init_valu1_fwrd, init_valu2_fwrd, init_valu1_bwrd, init_valu2_bwrd
end
end