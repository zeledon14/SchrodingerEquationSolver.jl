using AutomaticDocstrings

module PulayDensity
    """ a module to compute Pulay density mixing
    """ 

    mutable struct PulayData
        dns_mtrx::Vector{Vector{Float64}} #stored previous densities
        rsdl_mtrx::Vector{Vector{Float64}} #stored previous residuals
        R_mtrx::Matrix{Float64} #residual overlap matrix after dot product
        m::Int64 #maximum number of stored densities/residuals
        N::Int64 #number of grid points
        alpha::Float64 #mixing parameter
    end

    function init_pulay_data(N::Int64, m::Int64)::PulayData
        dns_mtrx= [zeros(Float64, N) for i in 1:m];
        rsdl_mtrx= [zeros(Float64, N) for i in 1:m];
        R_mtrx= zeros(Float64, m, m);
        return PulayData(dns_mtrx, rsdl_mtrx, R_mtrx, m, N, alpha);
    end


    function new_density_in(density_out::Vector{Float64, pula_data::PulayData, scl::Int64})::Vector{Float64}

        return new_density;

    end
end