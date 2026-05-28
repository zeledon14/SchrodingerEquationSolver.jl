module IntegralNumericalMethods

    """
    integrate_second_order_DE_RK4_PCABM5(grid::Vector{Float64}, 
        g::Vector{Float64}, f::Vector{Float64}, 
        init_valu1::Float64, init_valu2::Float64)::Vector{Float64}

        Solves differential equations of the form
            \frac{d^2y}{dx^2} = f(x)y + g(x) by transforming it into a system
            of equations.
            \frac{dy^0}{dx} = y^1(x) \\
            \frac{dy^1}{dx} = f(x)y^0(x) + g(x)
        to solve the system of differential equations, two methods are used Runge Kutta order 4
        (RK4) and prediction correction Adams-Moulton oder 5 (PCAM5). The first 4 values of 
        the solution uses the RK4 method and the rest uses the PCAM5.
            **Inputs:**
                -grid::Vector{Float64} values of the space grid where the functions are defined
                -g::Vector{Float64} a vector with the values g over the grid
                -f::Vector{Float64} a vector with the values f over the grid
                -init_valu1::Float64 first initial value of y at grid[1]
                -init_valu2::Float64 first initial value of y at grid[2]
                
            
            **Output:**
                -Vector{Float64} the function that solves the differential equation over the grid
"""
    #function integrate_second_order_DE_RK4_PCABM5_inital_cond_u
    function integrate_second_order_DE_RK4_PCABM5(grid::Vector{Float64}, 
        g::Vector{Float64}, f::Vector{Float64}, 
        init_valu1::Float64, init_valu2::Float64)::Vector{Float64}
        
        N=size(grid)[1];
        u=zeros(Float64, N);#solution to differential equation
        w=zeros(Float64, N);#first derivative of solution to differential equation
        u[1]= init_valu1;
        w[1]= (init_valu2 - init_valu1)/(grid[2]-grid[1]);

        for i in 1:4
        #for i in 1:(N-1)
            h= grid[i+1] - grid[i];
            u[i+1], w[i+1]= RK4(g[i:i+1],f[i:i+1],u[i], w[i],h);
        end
        #integration loop using prediction correction adams moulton degree 5
        for i in 6:N
            h= grid[i] - grid[i-1];
            u[i], w[i]= PCABM5(g[i-5:i],f[i-5:i],u[i-5:i-1], w[i-5:i-1],h);
        end
        
        return u
    end

    #function integrate_second_order_DE_RK4_PCABM5_inital_cond_u_w
    function integrate_second_order_DE_RK4_PCABM5_on_grid_i(grid_i::Vector{Float64}, 
        g::Vector{Float64}, f::Vector{Float64}, dr_di::Vector{Float64}, 
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
            
            u[i+1], w[i+1]= RK4_grid_i(g[i:i+1],f[i:i+1],dr_di[i:i+1],u[i], w[i],h);
        end
        #integration loop using prediction correction adams moulton degree 5
        for i in 6:N
            h=grid_i[i] -grid_i[i-1];
            #h=dr_di[i-1]*h;
            u[i], w[i]= PCABM5_grid_i(g[i-5:i],f[i-5:i], dr_di[i-5:i],
            u[i-5:i-1], w[i-5:i-1],h);
        end
        
        return u
    end

    #function integrate_second_order_DE_RK4_PCABM5_inital_cond_u_w
    function integrate_second_order_DE_RK4_PCABM5_direct_initial(grid::Vector{Float64}, 
        g::Vector{Float64}, f::Vector{Float64}, 
        u_in::Float64, du_in::Float64)::Vector{Float64}
        #direct initial because values of the 2 initial conditions 
        #are given directly. In contrast to the integrate_second_order_DE_RK4_PCABM5 
        #where the initial conditions are only direct for u and not for w
        N=size(grid)[1];
        u=zeros(Float64, N);#solution to differential equation
        w=zeros(Float64, N);#first derivative of solution to differential equation
        u[1]= u_in;
        w[1]= du_in;
        #println("function initial value", u[1])
        #println("function derivative initial value", w[1])
        #println("+++++++++++++++++++++++++++++++++++++")

        for i in 1:4
        #for i in 1:(N-1)
            h= grid[i+1] - grid[i];
            u[i+1], w[i+1]= RK4(g[i:i+1],f[i:i+1],u[i], w[i],h);
        end
        #integration loop using prediction correction adams moulton degree 5
        for i in 6:N
            h= grid[i] - grid[i-1];
            u[i], w[i]= PCABM5(g[i-5:i],f[i-5:i],u[i-5:i-1], w[i-5:i-1],h);
        end
        
        return u
    end

    function integrate_second_order_DE_Numerov(grid::Vector{Float64}, 
        g::Vector{Float64}, f::Vector{Float64}, 
        init_valu1::Float64, init_valu2::Float64)::Vector{Float64}
        
        N=size(grid)[1];
        u=zeros(Float64, N);#solution to differential equation

        u[1]= init_valu1;
        u[2]= init_valu2
        
        #for i in 1:4
        for i in 3:(N-2)
            h= grid[i+1] - grid[i];
            u[i],= Numerov(g[i-1],f[i-1],u[i-2:i-1],h);
        end
        
        return u
    end

"""
    RK4(g::Vector{Float64}, f::Vector{Float64},
         u::Float64, w::Float64, h::Float64)::Vector{Float64}

    Returns the values for y^0 and y^1 at the point x_{i+1} using
    Runge Kutta method of order 4 to solve equations of the form
    \frac{d^2y}{dx^2} = f(x)y + g(x) by transforming it into a system
    of equations.
    \frac{dy^0}{dx} = y^1(x) \\
    \frac{dy^1}{dx} = f(x)y^0(x) + g(x)

    **Inputs:**
        -g::Vector{Float64} a vector with the values [g_i, g_{i+1}]
        -f::Vector{Float64} a vector with the values [f_i, f_{i+1}]
        -u::Float64 value of the y^0(x_i)
        -w::Float64 value of the y^1(x_i)
        -h::Float64) value of the current delta x_{i+1} - x_{i}
    
    **Output:**
        -Vector{Float64} [y^0(x_i+1) y^1(x_i+1)] the values evaluated at x_{i+1}
"""
    function RK4(g::Vector{Float64}, f::Vector{Float64},
         u::Float64, w::Float64, h::Float64)::Tuple{Float64, Float64}
        #1 stands for the element i in the arrays
        #2 stands for the element i+1 in the arrays
        k01=h*(w)
        k11=h*(f[1]*u + g[1])

        k02=h*(w+0.5*k11)
        k12=h*(0.5*(f[1] + f[2])*(u+0.5*k01) +0.5*(g[1] + g[2]))

        k03=h*(w+0.5*k12)
        k13=h*(0.5*(f[1] + f[2])*(u+0.5*k02) +0.5*(g[1] + g[2]))

        k04=h*((w + k13) )
        k14=h*(f[2]*(u + k03) + g[2])

        up= u + (1.0/6.0)*(k01 + 2.0*k02 + 2.0*k03 + k04)
        wp= w + (1.0/6.0)*(k11 + 2.0*k12 + 2.0*k13 + k14)

        return up, wp
        
    end

    function RK4_grid_i(g::Vector{Float64}, f::Vector{Float64}, dr_di::Vector{Float64},
         u::Float64, w::Float64, h::Float64)::Tuple{Float64, Float64}
        #1 stands for the element i in the arrays
        #2 stands for the element i+1 in the arrays
        k01=h*dr_di[1]*(w)
        k11=h*dr_di[1]*(f[1]*u + g[1])

        k02=h*0.5*(dr_di[1] + dr_di[2])*(w+0.5*k11)
        k12=h*0.5*(dr_di[1] + dr_di[2])*(0.5*(f[1] + f[2])*(u+0.5*k01) +0.5*(g[1] + g[2]))

        k03=h*0.5*(dr_di[1] + dr_di[2])*(w+0.5*k12)
        k13=h*0.5*(dr_di[1] + dr_di[2])*(0.5*(f[1] + f[2])*(u+0.5*k02) +0.5*(g[1] + g[2]))

        k04=h*dr_di[2]*(w + k13) 
        k14=h*dr_di[2]*(f[2]*(u + k03) + g[2])

        up= u + (1.0/6.0)*(k01 + 2.0*k02 + 2.0*k03 + k04)
        wp= w + (1.0/6.0)*(k11 + 2.0*k12 + 2.0*k13 + k14)

        return up, wp
        
    end
"""
    PCAM4(g::Vector{Float64}, f::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Vector{Float64}
        Returns the values for y^0 and y^1 at the point x_{i} using
        predictor corrector Adams-Moulton order 4 to solve equations of the form
            \frac{d^2y}{dx^2} = f(x)y + g(x) by transforming it into a system
            of equations.
            \frac{dy^0}{dx} = y^1(x) \\
            \frac{dy^1}{dx} = f(x)y^0(x) + g(x)
        
            **Inputs:**
                -g::Vector{Float64} a vector with the values [g_{i-4}, g_{i-3}, g_{i-2}, g_{i-1}]
                -f::Vector{Float64} a vector with the values [f_{i-4}, f_{i-3}, f_{i-2}, f_{i-1}]
                -u::Vector{Float64} a vector with the values [y^0_{i-4}, y^0_{i-3}, y^0_{i-2}, y^0_{i-1}]
                -w::Vector{Float64} a vector with the values [y^1_{i-4}, y^1_{i-3}, y^1_{i-2}, y^1_{i-1}]
                -h::Float64) value of the current delta x_{i} - x_{i-1}
            
            **Output:**
                -Vector{Float64} [y^0(x_i) y^1(x_i)] the values evaluated at x_{i}
"""
    function PCAM4(g::Vector{Float64}, f::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Tuple{Float64, Float64}

        yp0= u[4] + (h/24.0)*(55.0*(w[4]) -59.0*(w[3]) 
                                +37.0*(w[2]) -9.0*(w[1]))
        yp1= w[4] + (h/24.0)*(55.0*(u[4]*f[4] + g[4]) -59.0*(u[3]*f[3] + g[3]) 
                                +37.0*(u[2]*f[2] + g[2]) -9.0*(u[1]*f[1] + g[1]))

        yc0= u[4] + (h/24.0)*(9.0*(yp1) +19.0*(w[4]) 
                                -5.0*(w[3]) +(w[2]))
        yc1= w[4] + (h/24.0)*(9.0*(yp0*f[5] + g[5]) +19.0*(u[4]*f[4] + g[4]) 
                                -5.0*(u[3]*f[3] + g[3]) +(u[2]*f[2] + g[2]))

        return yc0, yc1
    end

    function PCAM4_grid_i(g::Vector{Float64}, f::Vector{Float64}, dr_di::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Tuple{Float64, Float64}

        yp0= u[4] + (h/24.0)*(55.0*(dr_di[4]*w[4]) -59.0*(dr_di[3]*w[3]) 
                                +37.0*(dr_di[2]*w[2]) -9.0*(dr_di[1]*w[1]))
        yp1= w[4] + (h/24.0)*(55.0*dr_di[4]*(u[4]*f[4] + g[4]) -59.0*dr_di[3]*(u[3]*f[3] + g[3]) 
                                +37.0*dr_di[2]*(u[2]*f[2] + g[2]) -9.0*dr_di[1]*(u[1]*f[1] + g[1]))

        yc0= u[4] + (h/24.0)*(9.0*dr_di[5]*(yp1) +19.0*(dr_di[4]*w[4]) 
                                -5.0*(dr_di[3]*w[3]) +(dr_di[2]*w[2]))
        yc1= w[4] + (h/24.0)*(9.0*dr_di[5]*(yp0*f[5] + g[5]) +19.0*dr_di[4]*(u[4]*f[4] + g[4]) 
                                -5.0*dr_di[3]*(u[3]*f[3] + g[3]) +dr_di[2]*(u[2]*f[2] + g[2]))

        return yc0, yc1
    end
    """
    PCABM5(g::Vector{Float64}, f::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Vector{Float64}
        Returns the values for y^0 and y^1 at the point x_{i} using
            predictor corrector Adams-Moulton order 5 to solve equations of the form
                \frac{d^2y}{dx^2} = f(x)y + g(x) by transforming it into a system
                of equations.
                \frac{dy^0}{dx} = y^1(x) \\
                \frac{dy^1}{dx} = f(x)y^0(x) + g(x)
            
                **Inputs:**
                    -g::Vector{Float64} a vector with the values [g_{i-5}, g_{i-4}, g_{i-3}, g_{i-2}, g_{i-1}]
                    -f::Vector{Float64} a vector with the values [f_{i-5}, f_{i-4}, f_{i-3}, f_{i-2}, f_{i-1}]
                    -u::Vector{Float64} a vector with the values [y^0_{i-5}, y^0_{i-4}, y^0_{i-3}, y^0_{i-2}, y^0_{i-1}]
                    -w::Vector{Float64} a vector with the values [y^1_{i-5}, y^1_{i-4}, y^1_{i-3}, y^1_{i-2}, y^1_{i-1}]
                    -h::Float64) value of the current delta x_{i} - x_{i-1}
                
                **Output:**
                    -Vector{Float64} [y^0(x_i) y^1(x_i)] the values evaluated at x_{i}
"""
function PCABM5(g::Vector{Float64}, f::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Tuple{Float64, Float64}
        yp0= u[5] + (h/720.0)*(1901.0*(w[5]) 
                                -2774.0*(w[4]) +2616.0*(w[3]) 
                                -1274.0*(w[2]) +251.0*(w[1]))

        yp1= w[5] + (h/720.0)*(1901.0*(u[5]*f[5] + g[5]) 
                                -2774.0*(u[4]*f[4] + g[4]) +2616.0*(u[3]*f[3] + g[3]) 
                                -1274.0*(u[2]*f[2] + g[2]) +251.0*(u[1]*f[1] + g[1]))


        yc0= u[5] + (h/720.0)*(251.0*(yp1) 
                                +646.0*(w[5]) -264.0*(w[4]) 
                                +106.0*(w[3]) -19.0*(w[2]))

        yc1= w[5] + (h/720.0)*(251.0*(yp0*f[6] + g[6]) 
                                +646.0*(u[5]*f[5] + g[5]) -264.0*(u[4]*f[4] + g[4]) 
                                +106.0*(u[3]*f[3] + g[3]) -19.0*(u[2]*f[2] + g[2]))


        return yc0, yc1
    end

    function PCABM5_grid_i(g::Vector{Float64}, f::Vector{Float64}, dr_di::Vector{Float64},
        u::Vector{Float64}, w::Vector{Float64}, h::Float64)::Tuple{Float64, Float64}
        yp0= u[5] + (h/720.0)*(1901.0*(dr_di[5]*w[5]) 
                                -2774.0*(dr_di[4]*w[4]) +2616.0*(dr_di[3]*w[3]) 
                                -1274.0*(dr_di[2]*w[2]) +251.0*(dr_di[1]*w[1]))

        yp1= w[5] + (h/720.0)*(1901.0*(dr_di[5]*(u[5]*f[5] + g[5])) 
                                -2774.0*(dr_di[4]*(u[4]*f[4] + g[4])) +2616.0*(dr_di[3]*(u[3]*f[3] + g[3])) 
                                -1274.0*(dr_di[2]*(u[2]*f[2] + g[2])) +251.0*(dr_di[1]*(u[1]*f[1] + g[1])))


        yc0= u[5] + (h/720.0)*(251.0*(dr_di[6]*yp1) 
                                +646.0*(dr_di[5]*w[5]) -264.0*(dr_di[4]*w[4]) 
                                +106.0*(dr_di[3]*w[3]) -19.0*(dr_di[2]*w[2]))

        yc1= w[5] + (h/720.0)*(251.0*dr_di[6]*(yp0*f[6] + g[6]) 
                                +646.0*dr_di[5]*(u[5]*f[5] + g[5]) -264.0*dr_di[4]*(u[4]*f[4] + g[4]) 
                                +106.0*dr_di[3]*(u[3]*f[3] + g[3]) -19.0*dr_di[2]*(u[2]*f[2] + g[2]))


        return yc0, yc1
    end

    function Numerov(g::Float64, f::Float64,
        u::Vector{Float64}, h::Float64)::Float64
        yp= (u[2] + (h^2.0)*g/12.0)/(1.0 - (h^2.0)*f/12.0)
        u_out= 2.0*u[2] - u[1] + (h^2.0)*(f*yp + g)
       return u_out  
   end

end
