module Potentials

    """
    angular_potential(l::Int64,grid::Vector{Float64})::Vector{Float64}

    Produces the angular potential dependent on the angular momentum quantum number.
    **Inputs:**
        - l: angular momentum quantum number.
        - grid: points of the space where the functions are defined.
    **Output:**
        - angu_pote: angular potential defined on the grid points.

    """
    function angular_potential(l::Int64,grid::Vector{Float64})::Vector{Float64}
        angu_pote= l.*(l.+1.0)./(2.0.*grid.^2.0)
        return angu_pote
    end

    """
    coulomb_potential(Z::Int64 ,grid::Vector{Float64})::Vector{Float64}

    Produces the coulomb potential due to the nuclei acting on the electrons.
    **Inputs:**
        - Z: atomic number
        - grid: points of the space where the functions are defined.
    **Output:**
        - coul_pote: coulomb potential defined on the grid points.
    """
    function coulomb_potential(Z::Int64,grid::Vector{Float64})::Vector{Float64}
        coul_pote=-1.0.*Z./grid
        return coul_pote
    end

    """
    harmoic_oscilator_potential(w::Float64,m::Float64,grid::Vector{Float64})::Vector{Float64}

    produces a harmoic_oscilator_potential depending on the angular frequency, mass, and grid.
    **Inputs:**
        - w: angular frequency.
        - m: mass
        - grid: points of the space where the functions are defined.
    **Output:**
        - harm_pote: the harmonic oscilator potential defined on the grid points.
    """
    function harmoic_oscilator_potential(w::Float64,m::Float64,grid::Vector{Float64})::Vector{Float64}
        harm_pote= 0.5.*m.*(w^2).*grid.^2
        return harm_pote
    end
end