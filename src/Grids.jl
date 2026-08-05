
module Grids
    
    """
    exponential_grid(r_max::Float64, Z::Int64)::Vector{Float64}

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
    function exponential_grid(r_max::Float64, Z::Int64; 
                              b::Float64=0.002304)::Tuple{Vector{Float64}, Vector{Float64}}
        a=(4.34e-6)/Float64(Z)
        #b=0.002304
        #b=0.000504
        #b=0.001504
        
        N=((log((r_max/a)+1.0)/b))
        exp_grid= [a*(exp(b*i) -1.0) for i = 1.0:(N+1.0)]
        dx_di= a*b*exp.(b.*(1.0:(N+1.0)));
        return exp_grid, dx_di
    end


    mutable struct ExponentialGrid
        grid::Vector{Float64} #the exponenetial grid
        grid_i::Vector{Float64} #uniform number grid
        grid_sqrt::Vector{Float64} # exponenetial grid squared
        dx_di::Vector{Float64}#the dr in terms of i for integrals
        a::Float64
        b::Float64
        N::Int64
    end

    function ExponentialGrid(r_max::Float64,Z::Int64; 
        b::Float64=0.002304)::ExponentialGrid
        #b::Float64=0.002304;
        a::Float64=(4.34e-6)/Float64(Z);
        grid, dx_di= exponential_grid(r_max, Z, b=b);
        grid_i=[Float64(i) for (i,_) in enumerate(grid)];
        grid_sqrt= grid.^2;
        N=size(grid)[1];
        return ExponentialGrid(grid, grid_i, grid_sqrt,dx_di,a,b,N)
        
    end

    
    
    function simple_exponential_grid(r_min::Float64, r_max::Float64, N::Int64)::Vector{Float64}
        c::Float64=r_max/r_min
        grid= [r_min*c^(i/N) for i=0:(N-1)];
        return grid;
    end
    
    
    """
    uniform_grid(r_min::Float64, r_max::Float64, N::Int64)::Vector{Float64}

    Produces a unifor grid that starts at r_min and ends at r_max, and has 
    N points.
    **Inputs:**
        - r_min: initial r of the grid.
        - r_max: final r of the grid.
        - N: number of points in the grid.
    **Output:**
        - grid: a uniform grid, where the distance between succesive points in constant.

    """
    function uniform_grid(r_min::Float64, r_max::Float64, N::Int64)::Vector{Float64}
        delta= (r_max - r_min)/(N-1);
        grid= [r_min + i*delta for i=0:(N-1)];
        return grid;
    end

    mutable struct UniformGrid
        grid::Vector{Float64} #the grid
        grid_i::Vector{Float64} #uniform number grid
        grid_sqrt::Vector{Float64} # exponenetial grid squared
        dx_di::Vector{Float64} #derivateive dr/di
        #dx_di::Vector{Float64}#the dr in terms of i for integrals
        #a::Float64
        #b::Float64
        N::Int64
    end

    function UniformGrid(r_min::Float64, 
        r_max::Float64, N::Int64)::UniformGrid
        grid= uniform_grid(r_min, r_max, N);
        grid_i=[Float64(i) for (i,_) in enumerate(grid)];
        grid_sqrt= grid.^2;
        dx_di=fill(grid[2] - grid[1], N);
        #dx_di=(a*b).*(exp.(b.*grid_i));
        return UniformGrid(grid, grid_i, grid_sqrt,dx_di,N)
        
    end
end