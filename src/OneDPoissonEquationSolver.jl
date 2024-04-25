using AutomaticDocstrings

module OneDPoissonEquationSolver

using ..IntegralNumericalMethods
using ..MathUtils

raw"""
    solver(E::Float64,init_valu1_fwrd::Float64,
    init_valu2_fwrd::Float64, init_valu1_bwrd::Float64,
    init_valu2_bwrd::Float64, v_effe::Vector{Float64},
    grid::Vector{Float64})::Tuple{Vector{Float64},Float64}

    Returns the normalized wave function (u(r)) of one dimensional shcrodinger equation of the form
    \ frac{-1}{2} \ frac{d^2 u(r)}{dr^2} + V_{effe} u(r) = E u(r) 
    For numerical stability and accuracy is integrated forward and backward and merged at the first 
    classical turning point. The return solution is normalized.
    It also returns the merging value to help in the finding of the energy eigenvalue.

**Inputs:**
- `E::Float64`: Value to use as energy eigenvalue
- `init_valu1_fwrd::Float64`: Value of u(grid[1]) to use as first initial value in forward solution.
- `init_valu2_fwrd::Float64`: Value of u(grid[2]) to use as second initial value in forward solution.
- `init_valu1_bwrd::Float64`: Value of u(grid[end]) to use as first initial value in backward solution.
- `init_valu2_bwrd::Float64`: Value of u(grid[end -1]) to use as second initial value in backward solution.
- `v_effe::Vector{Float64}`: Vector with the values of the effective potential over the grid.
- `grid::Vector{Float64}`: Vector with the grid values.
"""

function solver(Z::Int64, density::Vector{Float64},
    grid::Vector{Float64})::Vector{Float64}

    g::Vector{Float64}=(-4.0*pi).*density.*grid;
    #g::Vector{Float64}=(-1).*density.*grid;
    f::Vector{Float64}= zeros(Float64, size(g)[1]);
    init_valu1_fwrd::Float64=Float64(grid[1])
    init_valu2_fwrd::Float64=Float64(grid[2])


    #do forward integration of poisson to get U_hartree
    U_hartree= IntegralNumericalMethods.integrate_second_order_DE_RK4_PCABM5(grid,g,f,
    init_valu1_fwrd,init_valu2_fwrd);
    #set boudary condtion over U_hartree
    a= (Z - U_hartree[end])/grid[end]
    U_hartree= U_hartree .+ a.*grid
    #transform into V_hartree
    V_hartree=U_hartree./grid#[U_hartree[i]/xi for (i,xi) in enumerate(grid)]
    #return U_hartree
    return V_hartree
end

end