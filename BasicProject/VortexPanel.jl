using Plots, Printf, LinearAlgebra, DelimitedFiles, FLOWFoil

using .AirfoilTools

#psuedocode
#=
    -Organize panels with coordinates for where each panel is starting and ending
        - I will probably need to give each panel a vector id with an x and y coordinate for each end of the panel
    -Need functions for
        - sin(θi)
        - cos(θi)
        - sin(θj)
        - cos(θj)
        - Bij
        - Aij matrix
        - Force of each panel computed from the tangent Velocity at each panel which comptues pressure coefficient and from that forces
        - Compute the Lift and drag coefficient by adding all the forces together vectorally to get the total lift and drag using the equation for lift and drag coefficient to get said coefficients.

=#
mutable struct Panels
        Panel_ID::Vector{Int64}
        Panel_start_points::Vector{Vector{Float64}}
        Panel_end_points::Vector{Vector{Float64}}
end

"""
    panel_setup(
    data_file;
    graph = false,
    graph_filename = ""
    )

Creates panels from points on an airfoil

# Arguments:
- `x::Array{Float}` : X coordinate's of the airfoil
- `z::Array{Float}` : z coordinate's of the airfoil

# Keyword Arguments:
- `graph::Boolean = false` : If true it will create a graph of the panels, if false it will not
- `graph_filename::String = ""` : If graph is true, graph_filename is a string containing the desired file location for the graph

# Returns:
- `Airfoil::Array` : Array with each panel ID, left, and right coordinates.
"""
function panel_setup(
    x,
    z;
    graph = false,
    graph_filename = "",
    )
    #this is a test to test out the syntax of structs with vectors
    #=
    test_panels = Panels([1,2],[[1.1,2.2], [3.1, 2.2]],[[1.1,2.2], [3.1, 2.2]])
    println(test_panels.Panel_start_points[2][2])
    =#
    Test1 = NACA4(3.0, 10.0, 12.0, false)
    x,z = naca4(Test1)
    startpoints = Vector{Vector{Float64}}
    endpoints = Vector{Vector{Float64}}
    ID = Vector{Int64}
    ID = [0]
    startpoints = [[0.0, 0.0]]
    endpoints = [[0.0, 0.0]]
    for i = 1:length(x) - 1
        push!(ID, i)
        startpoints = vcat(startpoints, [x[i],z[i]])
        endpoints = vcat(endpoints, [x[i + 1], z[i + 1]])
    end

    return 1 #test_panels
end

#Creates NACA coordinates using Airfoil AirfoilTools

#testing
#=
Test1 = NACA4(3.0, 10.0, 12.0, false)
x,z = naca4(Test1)
x = Vector{Vector{Float64}}
x = [[0.0, 0.0]]
x = vcat([[1.0, 2.0]])
println(x)
=#
x,z = 1,2
panel_setup(x,z)
println(" ")