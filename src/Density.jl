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

function adaptive_linear_mixing(density_in::Vector{Float64}, 
    density_out::Vector{Float64}; alpha_zero::Float64=0.5, beta::Float64=5.0)::Vector{Float64}

    """
    Performs adaptive linear mixing for density updates in SCF.

    Arguments:
    - density_in: Previous density (array)
    - density_out: New density obtained from solving KS equations (array)
    - alpha_zero: Initial mixing parameter (default 0.5)
    - beta: Scaling factor for adaptive step size control (default 5.0)

    Returns:
    - density_new: Updated density
    - alpha_k: Adaptive mixing parameter
    """

    # Compute residual norm
    residual_norm = (sum((density_out .- density_in).^2.0))^0.5 

    # Compute adaptive mixing parameter
    alpha_k = alpha_zero / (1 + beta * residual_norm)

    # Update density
    density_new = (1 - alpha_k) .* density_in .+ alpha_k .* density_out

    return density_new#, alpha_k
end



end