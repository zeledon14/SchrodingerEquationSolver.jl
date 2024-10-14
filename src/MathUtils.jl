using AutomaticDocstrings
module MathUtils

    """
    indices_of_zeros_finder(func::Vector{Float64})::Vector{Int64}

    Return the indices where the function func changes sign, the index 
    returned is the one before the change in sign. If func is of the form
    V_effe - E, the returned indices mark the classical turning points.

    **Input:**
        -func::Vector{Float64} the function to find the turning points
    
    **Output:**
        -indi::Vector{Int64} a vector with the indices of the turning points
            such that func[indi[i]] has a different sign than func[indi[i+1]]
    """
function indices_of_zeros_finder(func::Vector{Float64})::Vector{Int64}
    indi= Int64[]
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
    rescale!(solution1::Vector{Float64}, 
                   solution2::Vector{Float64}, 
                   turning_point::Int64)::Tuple{Vector{Float64},Vector{Float64}}
rescales the absolute biggest function to the smallest function
suth that solution1(turning_point) = solution2(turning_point). 

**Inputs:**
- `solution1::Vector{Float64}`: Vector with the values of the function1 to be scaled
- `solution2::Vector{Float64}`: Vector with the values of the function2 to be scaled
- `turning_point::Int64`: Point at wich v_effe - E = 0, there may be many of this points
                          the one use is lower.
**Output:**
- `solution1, solution2::Tuple{Vector{Float64},Vector{Float64}}`: Rescaled solutions.
"""
function  rescale!(solution1::Vector{Float64}, 
                   solution2::Vector{Float64}, 
                   turning_point::Int64)::Tuple{Vector{Float64},Vector{Float64}}
    A1=solution1[turning_point];
    A2=solution2[turning_point];
    if  abs(A1) > abs(A2)
        solution1= (A2/A1).*solution1
    else
        solution2= (A1/A2).*solution2
    end
    #We assume the solution 1 is always positive at solution1[1] this 
    #has to do with the odd states, that when changing the polarity made 
    #the search for eigenvalues harder
    if solution1[1] < 0.0
        solution1= (-1.0).*solution1
        solution2= (-1.0).*solution2
    end
    return solution1, solution2
end

"""
    three_point_derivative(func::Vector{Float64}, grid::Vector{Float64}, turning_point::Float64)

calculate the derivative of func over grid at the point turning_point,
the derivative uses a three point derivative.

**Inputs:**
- `func::Vector{Float64}`: the function for which the derivative is calculated
- `grid::Vector{Float64}`: the grid over which the function is defined
- `turning_point::Int64`: the point where the derivative is calculated
"""
function three_point_derivative(func::Vector{Float64},
                                grid::Vector{Float64},
                                turning_point::Int64)::Float64
    nume= func[turning_point+1] - func[turning_point-1]
    deno= (grid[turning_point+1] - grid[turning_point]) + (grid[turning_point] - grid[turning_point-1])
    temp= nume/deno
    return temp
end

"""
merge_solutions(forward::Vector{Float64}, backward::Vector{Float64},
                         grid::Vector{Float64}, turning_point::Int64)
                         ::Tuple{Vector{Float64},Float64}

merges the backward and forward integrated solutions into a merged one,
it also calculates the difference between the derivatives at the 
merging (turning) point, the merging happens at the turning_point.

**Inputs:**
- `forward::Vector{Float64}`: vextor with the forward integrated solution.
- `backward::Vector{Float64}`: vector with the backward integreated solution.
- `grid::Vector{Float64}`: the grid over which the function is defined.
- `turning_point::Float64`: turning point where the solutions are merged.

**Inputs:**
- `u_merged::Vector{Float64}`: vector with the merged function.
- `merge_value::Float64`: difference between derivatives on the turning (merging) point.
"""
function merge_solutions(forward::Vector{Float64}, backward::Vector{Float64},
    grid::Vector{Float64}, turning_point::Int64)::Tuple{Vector{Float64},Float64}
u_merged=zeros(Float64, size(forward)[1]);
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
    integral(func::Vector{Float64},grid::Vector{Float64})::Float64

integral(func::Vector{Float64}, grid::Vector{Float64})
Returns the value of the integral of func from grid[1] to grid[end]
the integration is perform using the trapezoidal rule.
I= sum 0.5*(func(x_{i+1})+func(x_{i}))*(x_{i+1} - x_{i})
**Inputs:**
- `func::Vector{Float64}`: vector with the values of the function to integrate
- `grid::Vector{Float64}`: grid where the function is defined
"""
function integral(func::Vector{Float64},grid::Vector{Float64})::Float64
    #temp=grid[2:end] .- grid[1:end-1];
    #temp1=func[2:end] .+ func[1:end-1];
    #I= temp.*temp1
    #I= 0.5*I
    I= 0.5*sum((grid[2:end] .- grid[1:end-1]).*(func[2:end] .+ func[1:end-1]))
    #I= 0.5*sum(I)
    return I
end

"""
    normalize!(func::Vector{Float64},grid::Vector{Float64})::Vector{Float64}

Normlizes the function func such that int func(x)^2 dx = 1
returns the normalized function func 

**Inputs:**
- `func::Vector{Float64}`: vector with the values of the function to normalize
- `grid::Vector{Float64}`: grid where the function is defined
"""
function normalize!(func::Vector{Float64},grid::Vector{Float64})::Vector{Float64}
    func_sqrt::Vector{Float64}= func.^2.0
    I= integral(func_sqrt, grid)
    I=I^(0.5)
    out= func./I
    return out
end

function normalize_v!(v::Vector{Float64},grid_struc::Any)::Vector{Float64}

    a::Float64=grid_struc.a;
    b::Float64=grid_struc.b;
    grid_i::Vector{Float64}=grid_struc.grid_i;

    func_sqrt::Vector{Float64}= (a*b.*exp.((2.0*b).*grid_i)).*(v.^2.0);

    I= integral(func_sqrt, grid_struc.grid_i);
    I=I^(0.5)
    out= v./I
    return out
end

"""
    error_difference(pred::Vector{Float64}, targ::Vector{Float64})

measures the error between the predicted calculated function, and the target
actual function.

**Inputs:**
- `pred::Vector{Float64}`: calculated function
- `targ::Vector{Float64}`: target function 
"""
function error_difference(pred::Vector{Float64},targ::Vector{Float64})::Float64
    temp::Vector{Float64}= ((pred .- targ).^2).^0.5;
    out::Float64=sum(temp)/length(temp)
    return out
end


function derivative(func::Vector{Float64},grid::Vector{Float64})::Vector{Float64}
    numerator= func[2:end] .- func[1:end-1];
    denomiator= grid[2:end] .- grid[1:end-1];
    out= numerator ./ denomiator
    return out
end

function is_continuous_enough(func::Vector{Float64}, 
                              indx::Int64, N_max::Int64)::Bool
    out::Bool=false
    if indx > 1 && (indx + 2) <= N_max
        if abs(func[indx+1]-func[indx]) < 
            0.51*(abs(func[indx+1]-func[indx+2]) + abs(func[indx]-func[indx-1]));
            out= true
            #the 0.51 gives some space for a continuos enough check, change for a better method
        end
    else
        throw(DomainError("indx to close to ends for continuoity conditions
        to be checked"));
    end
    return out
end

function energy_integral_exponential_grid(grid_stru::Any, density::Vector{Float64},
                                          vp::Vector{Float64})::Float64

    return integral((vp.*density.*(grid_stru.grid_sqrt).*(grid_stru.dr_i)), (grid_stru.grid_i));
    
end

end
