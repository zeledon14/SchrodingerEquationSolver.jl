using AutomaticDocstrings

AutomaticDocstrings.options[:min_args]=1
AutomaticDocstrings.options[:arg_types_in_desc]= true
AutomaticDocstrings.options[:args_header]="**Inputs:**"
AutomaticDocstrings.options[:arg_types_in_header]= true

#add @autodoc on top of the function
#On Julia  repl
#include("src/AutomaticDocstringsSetup.jl")
#autodoc("path_to_file")
