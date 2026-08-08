module IntegralNumericalMethods



    function integrate_second_order_DE_RK4_PCABM5_on_integer_grid(grid_i::Vector{Float64}, 
        g::Vector{Float64}, f::Vector{Float64}, dx_di::Vector{Float64}, 
        u_in::Float64, du_in::Float64)::Vector{Float64}
        #direct initial because values of the 2 initial conditions 
        #are given directly. In contrast to the integrate_second_order_DE_RK4_PCABM5 
        #where the initial conditions are only direct for u and not for w
        N=size(grid_i)[1];
        u=zeros(Float64, N);#solution to differential equation
        w=zeros(Float64, N);#first derivative of solution to differential equation
        u[1]= u_in;
        w[1]= du_in;
        #println("function initial value", u[1])
        #println("function derivative initial value", w[1])
        #println("+++++++++++++++++++++++++++++++++++++")

        for i in 1:4
        #for i in 1:(N-1)
            h=grid_i[i+1] -grid_i[i];
            
            u[i+1], w[i+1]= RK4_grid_i(g[i:i+1],f[i:i+1],dx_di[i:i+1],u[i], w[i],h);
        end
        #integration loop using prediction correction adams moulton degree 5
        for i in 6:N
            h=grid_i[i] -grid_i[i-1];
            #h=dx_di[i-1]*h;
            u[i], w[i]= PCABM5_grid_i(g[i-5:i],f[i-5:i], dx_di[i-5:i],
            u[i-5:i-1], w[i-5:i-1],h);
        end
        
        return u
    end

    function RK4_grid_i(g::Vector{Float64}, f::Vector{Float64}, dx_di::Vector{Float64},
         u::Float64, w::Float64, h::Float64)::Tuple{Float64, Float64}
        #1 stands for the element i in the arrays
        #2 stands for the element i+1 in the arrays
        k01=h*dx_di[1]*(w)
        k11=h*dx_di[1]*(f[1]*u + g[1])

        k02=h*0.5*(dx_di[1] + dx_di[2])*(w+0.5*k11)
        k12=h*0.5*(dx_di[1] + dx_di[2])*(0.5*(f[1] + f[2])*(u+0.5*k01) +0.5*(g[1] + g[2]))

        k03=h*0.5*(dx_di[1] + dx_di[2])*(w+0.5*k12)
        k13=h*0.5*(dx_di[1] + dx_di[2])*(0.5*(f[1] + f[2])*(u+0.5*k02) +0.5*(g[1] + g[2]))

        k04=h*dx_di[2]*(w + k13) 
        k14=h*dx_di[2]*(f[2]*(u + k03) + g[2])

        up= u + (1.0/6.0)*(k01 + 2.0*k02 + 2.0*k03 + k04)
        wp= w + (1.0/6.0)*(k11 + 2.0*k12 + 2.0*k13 + k14)

        return up, wp
        
    end

    function PCABM5_grid_i(g::Vector{Float64}, f::Vector{Float64}, dx_di::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Tuple{Float64, Float64}
        """Predictor corrector Adams-Bashforth-Moulton order 5 method for 
        solving second order differential equations"""
        yp0= u[5] + (h/720.0)*(1901.0*(dx_di[5]*w[5]) 
                                -2774.0*(dx_di[4]*w[4]) +2616.0*(dx_di[3]*w[3]) 
                                -1274.0*(dx_di[2]*w[2]) +251.0*(dx_di[1]*w[1]))

        yp1= w[5] + (h/720.0)*(1901.0*(dx_di[5]*(u[5]*f[5] + g[5])) 
                                -2774.0*(dx_di[4]*(u[4]*f[4] + g[4])) +2616.0*(dx_di[3]*(u[3]*f[3] + g[3])) 
                                -1274.0*(dx_di[2]*(u[2]*f[2] + g[2])) +251.0*(dx_di[1]*(u[1]*f[1] + g[1])))


        yc0= u[5] + (h/720.0)*(251.0*(dx_di[6]*yp1) 
                                +646.0*(dx_di[5]*w[5]) -264.0*(dx_di[4]*w[4]) 
                                +106.0*(dx_di[3]*w[3]) -19.0*(dx_di[2]*w[2]))

        yc1= w[5] + (h/720.0)*(251.0*dx_di[6]*(yp0*f[6] + g[6]) 
                                +646.0*dx_di[5]*(u[5]*f[5] + g[5]) -264.0*dx_di[4]*(u[4]*f[4] + g[4]) 
                                +106.0*dx_di[3]*(u[3]*f[3] + g[3]) -19.0*dx_di[2]*(u[2]*f[2] + g[2]))


        return yc0, yc1
    end


end
