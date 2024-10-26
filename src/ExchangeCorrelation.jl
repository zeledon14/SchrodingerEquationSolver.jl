using AutomaticDocstrings

module ExchangeCorrelation

    function potentials(density::Vector{Float64})::Tuple{Vector{Float64},
        Vector{Float64},Vector{Float64},Vector{Float64}}
        p=density
        A=0.0621814
        x0=-0.10498
        b=3.72744
        c=12.9352

        V_x= -1.0.*((3.0/pi).*p).^(1.0/3.0)
        E_x= (3.0/4.0).*V_x
        rs= (3.0./(4.0*pi.*p)).^(1.0/3.0)
        x= (rs).^0.5
        X_x= x.^2.0 .+ b.*x .+ c
        X_x0= x0.^2.0 .+ b.*x0 .+ c
        Q=(4.0*c-b^2.0)^0.5

        temp5= log.((x.^2.0)./(X_x))
        temp4= (2.0*b/Q).*atan.(Q./(2.0.*x.+b))
        temp3= (b*x0./X_x0)
        temp2= log.(((x.-x0).^2.0)./X_x)
        temp1= 2.0.*((b+2.0*x0)/(Q)).*(atan.(Q./(2.0.*x.+b)))
        temp6= (x .- x0) .- (b.*x0.*x)
        temp7= (x .- x0).*(X_x)

        E_c= replace_last_nans!((0.5*A).*(temp5 .+ temp4 .- temp3.*(temp2 .+ temp1)))
        V_c= E_c .- replace_last_nans!((A*c/6.0).*(temp6./temp7))
        
        return V_x, E_x, V_c, E_c
    end

    function replace_last_nans!(vec::Vector{Float64})::Vector{Float64}
        last_valid = findlast(!isnan, vec)  # Find the index of the last non-NaN element
        if last_valid === nothing
            return vec  # If no valid element is found, return the original vector
        end
        for i in (last_valid + 1):length(vec)
            vec[i] = vec[last_valid]
        end
        return vec
    end
    
end