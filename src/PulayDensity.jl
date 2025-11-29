using AutomaticDocstrings

module PulayDensity
    using ..MathUtils: integral
    """ a module to compute Pulay density mixing
    """ 

    mutable struct PulayData
        dns_mtrx::Matrix{Float64}#Vector{Vector{Float64}} #stored previous densities
        rsdl_mtrx::Matrix{Float64}#Vector{Vector{Float64}} #stored previous residuals
        R_mtrx::Matrix{Float64} #residual overlap matrix after dot product
        m::Int64 #maximum number of stored densities/residuals
        N::Int64 #number of grid points
        alpha::Float64 #mixing parameter
        dr_i::Vector{Float64}#the dr in terms of i for integrals
        grid_i::Vector{Float64} #uniform number grid
        grid_sqrt::Vector{Float64} # exponenetial grid squared
    end

    function init_pulay_data(N::Int64, m::Int64, alpha::Float64, dr_i, grid_i, grid_sqrt)::PulayData
        dns_mtrx= zeros(Float64, m, N);#[zeros(Float64, N) for i in 1:m];
        rsdl_mtrx= zeros(Float64, m, N);#[zeros(Float64, N) for i in 1:m];
        R_mtrx= zeros(Float64, m, m);
        return PulayData(dns_mtrx, rsdl_mtrx, R_mtrx, m, N, alpha, dr_i, grid_i, grid_sqrt);
    end

    function update_dns_rsdl_mtrx!(curr_dnst_in::Vector{Float64}, curr_dnst_out::Vector{Float64}, 
        pula_data::PulayData)
        # Update density matrix
        pula_data.dns_mtrx[2:end,:] = pula_data.dns_mtrx[1:end-1,:]
        pula_data.dns_mtrx[1,:] = curr_dnst_in
        # Update residual matrix
        pula_data.rsdl_mtrx[2:end, :] = pula_data.rsdl_mtrx[1:end-1, :]
        pula_data.rsdl_mtrx[1,:] = curr_dnst_out .- curr_dnst_in
    end

    function update_R_mtrx!(pula_data::PulayData, scl::Int64)   
        if scl == 1
            temp_prod= pula_data.rsdl_mtrx[1,:].* pula_data.rsdl_mtrx[1,:];
            pula_data.R_mtrx[1,1]= 4.0*pi*integral(temp_prod .* pula_data.grid_sqrt.*(pula_data.dr_i), pula_data.grid_i); 
        else
            L::Int64= min(scl, pula_data.m);
            #build temporal R temp_R
            temp_R::Vector{Float64}= zeros(Float64, L);
            for i in 1:L
                temp_prod= pula_data.rsdl_mtrx[1,:].* pula_data.rsdl_mtrx[i,:];
                temp_R[i]= 4.0*pi*integral(temp_prod .* pula_data.grid_sqrt.*(pula_data.dr_i), pula_data.grid_i); 
            end
            pula_data.R_mtrx[2:L,2:L] = pula_data.R_mtrx[1:L-1,1:L-1];
            pula_data.R_mtrx[1,1:L] = temp_R[1:L];
            pula_data.R_mtrx[1:L,1] = temp_R[1:L];
        end

    end

    function get_new_density_in(density_in::Vector{Float64}, density_out::Vector{Float64}, 
        pula_data::PulayData, scl::Int64)::Vector{Float64}
        L::Int64= min(scl, pula_data.m);
        update_dns_rsdl_mtrx!(density_in, density_out, 
        pula_data);
        update_R_mtrx!(pula_data, scl);
        #build B1 matrix
        B_1mtrx = ones(Float64, L+1, L+1);
        B_1mtrx[1:L,1:L] = pula_data.R_mtrx[1:L,1:L];
        B_1mtrx[L+1,L+1]= 0.0;
        println("B matrix:");
        println(B_1mtrx);
        zero_minus_one= zeros(Float64, L+1);
        zero_minus_one[end] = 1.0;
        B_1mtrx_inv= inv(B_1mtrx);
        println("B_1mtrx_inv:");
        println(B_1mtrx_inv);
        println("zero_minus_one")
        print(zero_minus_one);
        c_vector= B_1mtrx_inv * zero_minus_one;
        println("c_vector:");
        println(c_vector);
        #build new density
        new_density= zeros(Float64, pula_data.N);
        c_sum=0.0;
        for i in 1:L
            new_density .+= c_vector[i] .* (pula_data.dns_mtrx[i,:] .+ pula_data.alpha .* pula_data.rsdl_mtrx[i,:]);
            c_sum += c_vector[i];
        end
        println("c_sum:");
        println(c_sum);
        return new_density;

    end
end