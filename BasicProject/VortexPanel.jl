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
        Panel_mid_points::Vector{Vector{Float64}}
        Panel_length::Vector{Float64}
        theta::Vector{Float64}
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
- `panel_data::Panels` : Struct with each panel ID, left, right, and midpoint coordinates. Also includes length coordinates and angle coordinates relative to the x axis.
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

    #Initialize the vectors I will put in the Panels struct
    ID = [0]
    startpoints = [[0.0, 0.0]]
    endpoints = [[0.0, 0.0]]
    midpoints = [[0.0, 0.0]]
    l = [0.0]
    θ = [0.0]

    #add data to the vectors that will create the Panels struct
    for i = 1:length(x) - 1
        push!(ID, i)
        startpoints = vcat(startpoints, [[x[i],z[i]]])
        endpoints = vcat(endpoints, [[x[i + 1], z[i + 1]]])
        midpoints = vcat(midpoints, [[(x[i] + x[i + 1])/2 , (z[i] + z[i + 1])/2]])
        push!(l, sqrt((x[i + 1] - x[i])^2 + (z[i + 1] - z[i])^2))
        push!(θ, asin((z[i + 1] - z[i])/2))
    end
    #Get rid of the 0's that I put in to initialize the vectors I'm using
    splice!(ID, 1)
    splice!(startpoints, 1)
    splice!(endpoints, 1)
    splice!(midpoints, 1)
    splice!(l, 1)
    splice!(θ, 1)

    #create the struct
    panel_data = Panels(ID, startpoints, endpoints, midpoints, l, θ)

    #if graph is set to true then plot the airfoil
    if graph == true
        plot1 = plot(x,z)
        plot!(legend = false)
        xlims!(0.0, 1.0)
        ylims!(-0.5,0.5)
        savefig(plot1, graph_filename)
    end

    return panel_data
end

"""
    Beta_computer(
    Panel_data
    )

Computes Betaij values for an airfoil

# Arguments:
- `Panel_data::Panels` : Panels struct for the given airfoil

# Returns:
- `Betaij::Array` : Array with βij values, see Computational Aerodynamics by Andrew Ning equation 2.212
"""
function Beta_computer(
    Panel_data
)
    Betaij = Array{Float64, 2}(undef, length(Panel_data.Panel_ID), length(Panel_data.Panel_ID)) #create n x n Matrix

    #using equation 2.212 this for loop computes βij
    for i = 1:length(Panel_data.Panel_ID)
        for j = 1:length(Panel_data.Panel_ID)
            if i == j
                Betaij[i,j] = π*1
            else
                Betaij[i,j] = atan(
                        ((Panel_data.Panel_start_points[j][1] - Panel_data.Panel_mid_points[i][1])*(Panel_data.Panel_end_points[j][2] - Panel_data.Panel_mid_points[i][2])
                         - (Panel_data.Panel_start_points[j][2] - Panel_data.Panel_mid_points[i][2])*(Panel_data.Panel_end_points[j][1] - Panel_data.Panel_mid_points[i][1]))
                         ,
                            ((Panel_data.Panel_start_points[j][1] - Panel_data.Panel_mid_points[i][1])*(Panel_data.Panel_end_points[j][1] - Panel_data.Panel_mid_points[i][1])
                            + (Panel_data.Panel_start_points[j][2] - Panel_data.Panel_mid_points[i][2])*(Panel_data.Panel_end_points[j][2] - Panel_data.Panel_mid_points[i][2]))                       
                    )
            end
        end
    end
    return Betaij
end

"""
    Aij_computer(
    Panel_data,
    Betaij
    )

Computes Aij values from an airfoil

# Arguments:
- `Panel_data::Panels` : Panels struct for the given airfoil
- `Betaij::Array` : Array with βij values, see Computational Aerodynamics by Andrew Ning equation 2.212

# Returns:
- `Aij::Array` : Array with Aij values, see Computational Aerodynamic by Andrew Ning 2.223, 2.233, and 2.234
"""
function Aij_computer(
    Panel_data,
    Betaij
    )
    #initialize Aij Matrix to be an n + 1 x n + 1 matrix
    Aij = Matrix{Float64}(I, length(Panel_data.Panel_ID) + 1, length(Panel_data.Panel_ID) + 1)
    rij_1 = Matrix{Float64}(I, length(Panel_data.Panel_ID), length(Panel_data.Panel_ID))
    rij = Matrix{Float64}(I, length(Panel_data.Panel_ID), length(Panel_data.Panel_ID))
    #Computes n x n entries
    for i = 1:length(Panel_data.Panel_ID)
        for j = 1:length(Panel_data.Panel_ID)
            rij_1[i,j] = sqrt((Panel_data.Panel_end_points[j][1] - Panel_data.Panel_mid_points[i][1])^2 + (Panel_data.Panel_end_points[j][2] - Panel_data.Panel_mid_points[i][2])^2) #computes rij + 1
            rij[i,j] = sqrt((Panel_data.Panel_start_points[j][1] - Panel_data.Panel_mid_points[i][1])^2 + (Panel_data.Panel_start_points[j][2] - Panel_data.Panel_mid_points[i][2])^2) #computes rij
            Aij[i, j] = (
                log(ℯ, rij_1[i,j] / rij[i,j])*sin(Panel_data.theta[i] - Panel_data.theta[j]) + Betaij[i,j]*cos(Panel_data.theta[i] - Panel_data.theta[j]) #computes Aij
            )
        end
    end

    #computes [1 through n rows, n + 1 column] for Aij see equation 2.223
    for i = 1:length(Panel_data.Panel_ID)
        for j = 1:length(Panel_data.Panel_ID)
            Aij[i, length(Panel_data.Panel_ID) +  1] = Aij[i, length(Panel_data.Panel_ID) +  1] + log(ℯ, (rij_1[i,j] / rij[i,j]))*cos(Panel_data.theta[i] - Panel_data.theta[j]) - Betaij[i,j]*sin(Panel_data.theta[i] - Panel_data.theta[j])
        end
    end

    #computes [n + 1 row, 1 through n columns]
    for j = 1:length(Panel_data.Panel_ID)
        Aij[length(Panel_data.Panel_ID) + 1, j] = (
            Betaij[1,j]*sin(Panel_data.theta[1] - Panel_data.theta[j]) - log(ℯ, (rij_1[1,j] / rij[1,j]))*cos(Panel_data.theta[1] - Panel_data.theta[j])
            + Betaij[length(Panel_data.Panel_ID),j]*sin(Panel_data.theta[length(Panel_data.Panel_ID)] - Panel_data.theta[j]) - log(ℯ, (rij_1[length(Panel_data.Panel_ID),j] / rij[length(Panel_data.Panel_ID),j]))*cos(Panel_data.theta[length(Panel_data.Panel_ID)] - Panel_data.theta[j])
            )
    end

    #computes [n + 1 row, n + 1 column]
    for j = 1:length(Panel_data.Panel_ID)
        Aij[length(Panel_data.Panel_ID) + 1, length(Panel_data.Panel_ID) + 1] = (
            Aij[length(Panel_data.Panel_ID) + 1, length(Panel_data.Panel_ID) + 1]
            + Betaij[1,j]*cos(Panel_data.theta[1] - Panel_data.theta[j])
            + log(ℯ, (rij_1[1,j] / rij[1,j]))*sin(Panel_data.theta[1] - Panel_data.theta[j])
            ) + (
                Betaij[length(Panel_data.Panel_ID),j]*cos(Panel_data.theta[length(Panel_data.Panel_ID)] - Panel_data.theta[j])
                +  log(ℯ, (rij_1[length(Panel_data.Panel_ID),j] / rij[length(Panel_data.Panel_ID),j]))*sin(Panel_data.theta[length(Panel_data.Panel_ID)] - Panel_data.theta[j])
            )
    end

    return Aij
end
#Creates NACA coordinates using Airfoil AirfoilTools
Test1 = NACA4(2.0, 4.0, 12.0, false)
x,z = naca4(Test1)
xvalues = [0.0]

#call panel setup function
test_panels = panel_setup(x,z, graph = true, graph_filename = "BasicProject\\TestGraph.png")
Betaij = Beta_computer(test_panels)
Aij = Aij_computer(test_panels, Betaij)
println(" ")