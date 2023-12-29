
module Grids
    
    """
    exponential_grid(r_max::Float32, Z::Int32)::Vector{Float32}

    Exponential grid as defined in J.P. Desclaux, Comp. Phys. Comm. 1, 216 (1969).
    Reference to the data to build the grid can be access in 
    https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-3
    **Inputs:**
        - r_max: maximal radius of the grid.
        - Z: atomic number of th element for which the grid is constructed.
    **Output:**
        - exp_grid: this grid does no includes 0 exp_grid[0] depends of Z and
            0 is excluded to avoid infinits. The grid goes up to r_max. The number
            of elements in the grid is control internaly.
    
    """
    function exponential_grid(r_max::Float32, Z::Int32)::Vector{Float32}
        a=(4.34*10.0^(-6.0))/Z
        b=0.002304
        
        N=Integer(ceil(log((r_max/a)+1.0)/b))

        exp_grid= [a*(exp(b*i) -1.0) for i = 1:(N+1)]
        return exp_grid 
    end


    """
    uniform_grid(r_min::Float32, r_max::Float32, N::Int32)::Vector{Float32}

    Produces a unifor grid that starts at r_min and ends at r_max, and has 
    N points.
    **Inputs:**
        - r_min: initial r of the grid.
        - r_max: final r of the grid.
        - N: number of points in the grid.
    **Output:**
        - grid: a uniform grid, where the distance between succesive points in constant.

    """
    function uniform_grid(r_min::Float32, r_max::Float32, N::Int32)::Vector{Float32}
        delta= (r_max - r_min)/(N-1)
        grid= [r_min + i*delta for i=0:(N-1)]
        return grid
    end
end