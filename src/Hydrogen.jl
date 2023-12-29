module Hydrogen

    """
    u_s1_hydrogen(grid::Vector{Float32})::Vector{Float32}

    radial u function that is the solution of the radial shcrodinger equation  
    in the s1 state with eigenvalue 0.5
    
    **Inputs:**
        - grid: The vector of points where the function is calcuated.

    **Output:**
        - u_s1: The vector with the values of u_s1 over the grid.
    """
    function u_s1_hydrogen(grid::Vector{Float32})::Vector{Float32}
        u=[xi*exp(-1.0*xi) for xi in grid]
        return u 
    end

    function U_hartree(grid)
        out=[(1.0 -1.0*(xi+1.0)*exp(-2.0*xi))  for (i,xi) in enumerate(grid)]
    end
end