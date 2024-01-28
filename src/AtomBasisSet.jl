using AutomaticDocstrings


module AtomBasisSet
using JSON3

    
mutable struct orbital
    n::Int64;
    l::Int64;
    name::String;
    E::Float64;
    u:: Vector{Float64};
end

mutable struct atom_basis_set
    grid::Vector{Float64}
    orbitals::Array{orbital}
end



function init_atom_basis_set(Z::Int64, grid::Vector{Float64})::atom_basis_set
    file= open("../data/atomic_numbers.json");
    atomic_numbers= JSON3.read(file);
    N::Int64=sizeof(grid)[1];
    electron_capacity::Int64=0;
    numb_orbitals::Int64=0
    #finding the number of orbitals needed given Z
    for i_orbi in atomic_numbers
        if electron_capacity + i_orbi.max_number_of_electrons <= Z
            electron_capacity += i_orbi.max_number_of_electrons;
            numb_orbitals+=1;
        end  
    end
    #init orbitals
    #initl with the u wave function as zeroes
    orbitals=[orbital(0,0,"",0,zeros(N)) for _ in (1:numb_orbitals)]
    #fill orbitals with information
    for i in (1:numb_orbitals) 
        orbitals[i].n=atomic_numbers[i].n
        orbitals[i].l=atomic_numbers[i].l
        orbitals[i].name=atomic_numbers[i].orbital_name
        orbitals[i].E= -0.5*(Z^2)/(atomic_numbers[i].n^2)
    end

    return atom_basis_set(grid,orbitals)
    
end
end