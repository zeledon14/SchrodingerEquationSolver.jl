using AutomaticDocstrings

module Density

function  calculate_density(basis_set::Any)::Vector{Float64}
    density::Vector{Float64}= zero(basis_set.grid)
    grid::Vector{Float64}= basis_set.grid
    for i_orbi in basis_set.orbitals
        c= i_orbi.occu/(4.0*pi)
        #c= i_orbi.occu
        density .= density .+ c.*(i_orbi.u.^2.0)./(grid.^2.0)
    end
    return density
end    

function linear_mixing(density_in::Vector{Float64}, 
    density_out::Vector{Float64}; alpha::Float64=0.3)::Vector{Float64}
    temp1= alpha.*density_out
    temp2= (1.0 - alpha).*density_in
    return temp1 .+ temp2
end


end