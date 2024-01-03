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

function harmoic_oscillator(grid::Vector{Float64},  E::Float64=-0.5, 
    l::Int64=0)::Tuple{Float64,Float64,Float64,Float64}
    init_valu1_fwrd::Float64=exp(-0.5*abs(grid[1])^2);
    init_valu2_fwrd::Float64=exp(-0.5*abs(grid[2])^2);
    init_valu1_bwrd::Float64=exp(-0.5*abs(grid[end])^2);
    init_valu2_bwrd::Float64=exp(-0.5*abs(grid[end-1])^2);
    return init_valu1_fwrd, init_valu2_fwrd, init_valu1_bwrd, init_valu2_bwrd
end
end