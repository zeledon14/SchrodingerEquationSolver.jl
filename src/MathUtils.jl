using AutomaticDocstrings
module MathUtils

    """
    indices_of_zeros_finder(func::Vector{Float32})::Vector{Int32}

    Return the indices where the function func changes sign, the index 
    returned is the one before the change in sign. If func is of the form
    V_effe - E, the returned indices mark the classical turning points.

    **Input:**
        -func::Vector{Float32} the function to find the turning points
    
    **Output:**
        -indi::Vector{Int32} a vector with the indices of the turning points
            such that func[indi[i]] has a different sign than func[indi[i+1]]
    """
function indices_of_zeros_finder(func::Vector{Float32})::Vector{Int32}
    indi= Int32[]
    sign_befo= Integer(sign(func[1]))
    for (i,f_i) in enumerate(func[2:end])
        if sign_befo != Integer(sign(f_i))
            append!(indi, i)
            sign_befo = Integer(sign(f_i))
        end
    end
    return indi
        
end

"""
    rescale!(solution1::Vector{Float32}, 
                   solution2::Vector{Float32}, 
                   turning_point::Int32)::Tuple{Vector{Float32},Vector{Float32}}
rescales the absolute biggest function to the smallest function
suth that solution1(turning_point) = solution2(turning_point). 

**Inputs:**
- `solution1::Vector{Float32}`: Vector with the values of the function1 to be scaled
- `solution2::Vector{Float32}`: Vector with the values of the function2 to be scaled
- `turning_point::Int32`: Point at wich v_effe - E = 0, there may be many of this points
                          the one use is lower.
**Output:**
- `solution1, solution2::Tuple{Vector{Float32},Vector{Float32}}`: Rescaled solutions.
"""
function  rescale!(solution1::Vector{Float32}, 
                   solution2::Vector{Float32}, 
                   turning_point::Int32)::Tuple{Vector{Float32},Vector{Float32}}
    A1=solution1[turning_point];
    A2=solution2[turning_point];
    if  abs(A1) > abs(A2)
        solution1= (A2/A1).*solution1
    else
        solution2= (A1/A2).*solution2
    end
    return solution1, solution2
end

"""
    three_point_derivative(func::Vector{Float32}, grid::Vector{Float32}, turning_point::Float32)

calculate the derivative of func over grid at the point turning_point,
the derivative uses a three point derivative.

**Inputs:**
- `func::Vector{Float32}`: the function for which the derivative is calculated
- `grid::Vector{Float32}`: the grid over which the function is defined
- `turning_point::Int32`: the point where the derivative is calculated
"""
function three_point_derivative(func::Vector{Float32},
                                grid::Vector{Float32},
                                turning_point::Int32)::Float32
    nume= func[turning_point+1] - func[turning_point-1]
    deno= (grid[turning_point+1] - grid[turning_point]) + (grid[turning_point] - grid[turning_point-1])
    temp= nume/deno
    return temp
end

"""
merge_solutions(forward::Vector{Float32}, backward::Vector{Float32},
                         grid::Vector{Float32}, turning_point::Int32)
                         ::Tuple{Vector{Float32},Float32}

merges the backward and forward integrated solutions into a merged one,
it also calculates the difference between the derivatives at the 
merging (turning) point, the merging happens at the turning_point.

**Inputs:**
- `forward::Vector{Float32}`: vextor with the forward integrated solution.
- `backward::Vector{Float32}`: vector with the backward integreated solution.
- `grid::Vector{Float32}`: the grid over which the function is defined.
- `turning_point::Float32`: turning point where the solutions are merged.

**Inputs:**
- `u_merged::Vector{Float32}`: vector with the merged function.
- `merge_value::Float32`: difference between derivatives on the turning (merging) point.
"""
function merge_solutions(forward::Vector{Float32}, backward::Vector{Float32},
                         grid::Vector{Float32}, turning_point::Int32)::Tuple{Vector{Float32},Float32}
    u_merged=zeros(Float32, size(forward)[1]);
    fwrd_drvt= three_point_derivative(forward, grid, turning_point);
    bwrd_drvt= three_point_derivative(backward, grid, turning_point);
    merge_value= fwrd_drvt -  bwrd_drvt;
    u_merged[1:turning_point-2]= forward[1:turning_point-2];
    u_merged[turning_point-1:turning_point+1]= (0.5).*(forward[turning_point-1:turning_point+1] 
                                                       .+ backward[turning_point-1:turning_point+1]);
    u_merged[turning_point+2:end]= backward[turning_point+2:end];
    return u_merged, merge_value
end

"""
    integral(func::Vector{Float32},grid::Vector{Float32})::Float32

integral(func::Vector{Float32}, grid::Vector{Float32})
Returns the value of the integral of func from grid[1] to grid[end]
the integration is perform using the trapezoidal rule.
I= sum 0.5*(func(x_{i+1})+func(x_{i}))*(x_{i+1} - x_{i})
**Inputs:**
- `func::Vector{Float32}`: vector with the values of the function to integrate
- `grid::Vector{Float32}`: grid where the function is defined
"""
function integral(func::Vector{Float32},grid::Vector{Float32})::Float32
    I= 0.5.*(grid[2:end] .- grid[1:end-1]).*(func[2:end] .+ func[1:end-1])
    I= sum(I)
    return I
end

"""
    normalize!(func::Vector{Float32},grid::Vector{Float32})::Vector{Float32}

Normlizes the function func such that int func(x)^2 dx = 1
returns the normalized function func 

**Inputs:**
- `func::Vector{Float32}`: vector with the values of the function to normalize
- `grid::Vector{Float32}`: grid where the function is defined
"""
function normalize!(func::Vector{Float32},grid::Vector{Float32})::Vector{Float32}
    func_sqrt::Vector{Float32}= func.^2.0
    I= integral(func_sqrt, grid)
    I=I^(0.5)
    out= func./I
    return out
end

"""
    error_difference(pred::Vector{Float32}, targ::Vector{Float32})

measures the error between the predicted calculated function, and the target
actual function.

**Inputs:**
- `pred::Vector{Float32}`: calculated function
- `targ::Vector{Float32}`: target function 
"""
function error_difference(pred::Vector{Float32},targ::Vector{Float32})::Float32
    temp::Vector{Float32}= ((pred .- targ).^2).^0.5;
    out::Float32=sum(temp)/length(temp)
    return out
end

end
