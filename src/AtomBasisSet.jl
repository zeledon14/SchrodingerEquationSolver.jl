using AutomaticDocstrings


module AtomBasisSet
using JSON3

"""
   orbital


   Structure to hold information about a given radial orbital


# Fields:
- `n::Int64`: principal quantum number.
- `l::Int64`: orbital quantum number.
- `name::String`: orbital name.
- `E::Float64`: energy eigenvalue of the orbital .
- `occu::Float64`: number of electrons in the radial orbital.
- `u::Vector{Float64}`: related to the solution of the radial Schrodinger
                       equation u = r (Psi), where Psi is the solution
                       of the Schrodinger equation.
"""
mutable struct orbital
    n::Int64;
    l::Int64;
    name::String;
    E::Float64;
    occu::Float64;
    u:: Vector{Float64};
end

function orbital_to_dict(in_orbital::orbital)
    out= Dict(fieldnames(orbital) .=> getfield.(Ref(in_orbital), 
                                            fieldnames(orbital)));
    return out
end
"""
   atom_basis_set


   Structure to hold a list of orbitals and grid making a basis set.


# Fields:
- `grid::Vector{Float64}`: A list with the space points.
- `orbitals::Array{orbital}`: A list with the orbitals in the basis set.
"""
mutable struct atom_basis_set
    Name::String
    Energy::Float64
    E_hartree::Float64
    E_xc::Float64
    V_xc::Float64
    grid::Vector{Float64}
    orbitals::Array{orbital}
end

"""
   init_atom_basis_set(Z::Int64, grid::Vector{Float64})
**Inputs:**
- `Z::Int64`: Atomic number of the element.
- `grid::Vector{Float64}`: A list with the space points.
"""
function init_atom_basis_set(Z::Int64, grid::Vector{Float64})::atom_basis_set
    file= open(joinpath(dirname(@__FILE__),"../data/atomic_numbers.json"));
    atomic_numbers= JSON3.read(file);
    file1= open(joinpath(dirname(@__FILE__),"../data/atom_number_name.json"));
    Name= JSON3.read(file1)[Int64(Z)];
    N::Int64=size(grid)[1];
    #electron_capacity::Int64=0;
    numb_orbitals::Int64=0
    remaining_electrons::Int64= Int64(Z)
    #finding number of orbitals needed
    for i_orbi in atomic_numbers
        if remaining_electrons > 0
            remaining_electrons -= i_orbi.max_number_of_electrons;
            numb_orbitals+=1;
        end  
    end
    #init orbitals
    #initl with the u wave function as zeroes
    orbitals=[orbital(0,0,"",0,0,zeros(N)) for _ in (1:numb_orbitals)]
    remaining_electrons= Int64(Z)
    #fill orbitals with information
    for i in (1:numb_orbitals) 
        orbitals[i].n=atomic_numbers[i].n
        orbitals[i].l=atomic_numbers[i].l
        orbitals[i].name=atomic_numbers[i].orbital_name
        orbitals[i].E= -0.5*(Z^2)/(atomic_numbers[i].n^2)
        if remaining_electrons >= atomic_numbers[i].max_number_of_electrons
            orbitals[i].occu= atomic_numbers[i].max_number_of_electrons
            remaining_electrons -= atomic_numbers[i].max_number_of_electrons
        else
            orbitals[i].occu=remaining_electrons
        end

    end
    #Initial values
    Energy::Float64=0.0
    E_hartree::Float64=0.0
    E_xc::Float64=0.0
    V_xc::Float64=0.0
    return atom_basis_set(Name,Energy,E_hartree,
                            E_xc,V_xc,
                            grid,orbitals)
    
end

function save_basis_set(in_basis::atom_basis_set, save_path::String)    
    orbitals= [orbital_to_dict(i_orbi) for i_orbi in in_basis.orbitals];
    out= JSON3.write(Dict{String, Any}("grid"=>in_basis.grid,
                             "Name"=>in_basis.Name,
                             "Energy"=>in_basis.Energy,
                             "E_hartree"=>in_basis.E_hartree,
                             "E_xc"=>in_basis.E_xc,
                             "V_xc"=>in_basis.grid,
                             "orbitals"=>orbitals));
    open(save_path,"w") do f 
        write(f, out) 
    end 
end
function load_basis_set_from_json(path::String)
    return JSON3.read(path, atom_basis_set);
end
end