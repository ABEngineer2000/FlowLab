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
        sin_theta::Vector{Float64}
        cos_theta::Vector{Float64}
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
    sin_t = [0.0]
    cos_t = [0.0]

    #add data to the vectors that will create the Panels struct
    for i = 1:length(x) - 1
        push!(ID, i)
        startpoints = vcat(startpoints, [[x[i],z[i]]])
        endpoints = vcat(endpoints, [[x[i + 1], z[i + 1]]])
        midpoints = vcat(midpoints, [[(x[i] + x[i + 1])/2 , (z[i] + z[i + 1])/2]])
        push!(l, sqrt((x[i + 1] - x[i])^2 + (z[i + 1] - z[i])^2))
        push!(θ, atan(z[i + 1] - z[i], x[i + 1] - x[i]))
        push!(sin_t, (z[i + 1] - z[i]) / l[i + 1])
        push!(cos_t, (x[i + 1] - x[i]) / l[i + 1])
    end
    #Get rid of the 0's that I put in to initialize the vectors I'm using
    splice!(ID, 1)
    splice!(startpoints, 1)
    splice!(endpoints, 1)
    splice!(midpoints, 1)
    splice!(l, 1)
    splice!(θ, 1)
    splice!(sin_t, 1)
    splice!(cos_t, 1)

    #create the struct
    panel_data = Panels(ID, startpoints, endpoints, midpoints, l, θ, sin_t, cos_t)

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
- `sinij::Array` : sinij Array - see thetaij
- `cosij::Array` : cosij Array - see thetaij

# Returns:
- `Aij::Matrix` : Array with Aij values, see Computational Aerodynamics by Andrew Ning 2.223, 2.233, and 2.234
"""
function Aij_computer(
    Panel_data,
    Betaij,
    sinij,
    cosij
    )
    n = length(Panel_data.Panel_ID)
    #initialize Aij Matrix to be an n + 1 x n + 1 matrix
    Aij = Matrix{Float64}(I, n + 1, n + 1)
    rij_1 = Matrix{Float64}(I, n, n)
    rij = Matrix{Float64}(I, n, n)
    #Computes n x n entries
    for i = 1:n
        for j = 1:n
            rij_1[i,j] = sqrt((Panel_data.Panel_end_points[j][1] - Panel_data.Panel_mid_points[i][1])^2 + (Panel_data.Panel_end_points[j][2] - Panel_data.Panel_mid_points[i][2])^2) #computes rij + 1
            rij[i,j] = sqrt((Panel_data.Panel_start_points[j][1] - Panel_data.Panel_mid_points[i][1])^2 + (Panel_data.Panel_start_points[j][2] - Panel_data.Panel_mid_points[i][2])^2) #computes rij
            Aij[i, j] = (
                log(ℯ, rij_1[i,j] / rij[i,j])*sinij[i,j] + Betaij[i,j]*cosij[i,j]
                ) #computes Aij
            
        end
    end

    #computes [1 through n rows, n + 1 column] for Aij see equation 2.223
    for i = 1:n
        for j = 1:n
            Aij[i, n + 1] = Aij[i, n + 1] + log(ℯ, (rij_1[i,j] / rij[i,j]))*cosij[i,j] - Betaij[i,j]*sinij[i,j]
        end
    end

    #computes [n + 1 row, 1 through n columns]
    for j = 1:n
        Aij[n + 1, j] = (
            Betaij[1,j]*sinij[1,j] - log(ℯ, (rij_1[1,j] / rij[1,j]))*cosij[1,j]
            + Betaij[n,j]*sinij[n,j] - log(ℯ, (rij_1[n,j] / rij[n,j]))*cosij[n,j]
            )
    end

    #computes [n + 1 row, n + 1 column]
    for j = 1:n
        Aij[n + 1, n + 1] = (
            Aij[n + 1, n + 1]
            + Betaij[1,j]*cosij[1,j]
            + log(ℯ, (rij_1[1,j] / rij[1,j]))*sinij[1,j]
            ) + (
                Betaij[n,j]*cosij[n,j]
                +  log(ℯ, (rij_1[n,j] / rij[n,j]))*sinij[n,j]
            )
    end

    return Aij, rij, rij_1
end


"""
    B_computer(
    Panel_data,
    α,
    Vinf
    )

Computes B valeus from an airfoil see Computational Aerodynamics by Andrew Ning equations 2.233 and 2.232

# Arguments:
- `Panel_data::Panels` : Panels struct for the given airfoil
- `α::Float` : User specified angle of attack in degrees
- `Vinf::Float` : User specified free streamm velocity relative to the airfoil chord length

# Returns:
- `B::Matrix` : Matrix with B values
"""
function B_computer(
    Panel_data,
    α,
    Vinf,
    ) 
    α = α*π / 180 #convert α to radians
    n = length(Panel_data.Panel_ID)
    B = Matrix{Float64}(undef, n + 1, 1)
    
    #compute 1 through n rows
    for i = 1:n
        B[i,1] = 2*π*Vinf*(Panel_data.sin_theta[i]*cos(α) - Panel_data.cos_theta[i]*sin(α))      #sin(Panel_data.theta[i] - α)
    end

    #compute n +1 row
    B[n + 1,1] = -2*π*Vinf*(Panel_data.cos_theta[1]*cos(α) + Panel_data.sin_theta[1]*sin(α) + Panel_data.cos_theta[n]*cos(α) + Panel_data.sin_theta[n]*sin(α))            #(cos(Panel_data.theta[1] - α) + cos(Panel_data.theta[n]  - α))
    
    return B
end

"""
    solve_system(
    Aij,
    B
    )

Solves system of equations for the Hess-Smith Panel Method - see Computational Aerodynamics by Andrew Ning equation 2.234

# Arguments:
- `Aij::Matrix` : Aij Matrix - see Aij_computer
- `B::Matrix` : B values - see B_computer

# Returns:
- `solution::Matrix` : n + 1 x 1 Matrix where 1 through n rows are the source strengths and the n + 1 row is the vortex strength
"""
function solve_system(
    Aij,
    B
)
    n = length(B)
    solution = Matrix{Float64}(undef, n, 1)
    solution = inv(Aij)*B
    return solution
end

"""
    Tangent_Velocity_computer(
    Panel_data,
    Aij,
    βij,
    q_λ_vector,
    vinf,
    α
    )

Computes the tangent velocity at the control point of each of the panels

# Arguments:
- `Panel_data::Panels` : Panel struct containing panel data
- `rij::Matrix` : rij - see Aij_computer
- `rij_1::Matrix` : rij_1 - see Aij_computer
- `βij::Matrix` : βij Matrix - see Beta_computer
- `sinij::Array` : sinij Array - see thetaij
- `cosij::Array` : cosij Array - see thetaij
- `q_λ_vector::Matrix` : Solution vector - see solve_system
- `Vinf::Float` : Free stream velocity relative to the airfoil chord length
- `α::Float` : Angle of attack in degrees

# Returns:
- `Vti::Vector` : n length Vector which holds the value of the tangent velocity at the control point for each panel
"""
function Tangent_velocity_computer(
    Panel_data,
    rij,
    rij_1,
    βij,
    sinij,
    cosij,
    q_λ_vector,
    Vinf,
    α
    )
    α = α*π/180 #convert α to radians
    n = length(Panel_data.Panel_ID)
    Vti = Vector{Float64}(undef, n)

    #compute's tangent velocity at each of the panels see equation 2.237 in Computational Aerodynamics by Andrew Ning
    for i = 1:n
        for j = 1:n
            Vti[i] = (Vti[i] +
            (1 / (2*π))*q_λ_vector[j]*(
                βij[i,j]*sinij[i,j] - log(ℯ, (rij_1[i,j] / rij[i,j]))*cosij[i,j]
            )
            + (q_λ_vector[n + 1] / (2*π))*(
                βij[i,j]*cosij[i,j] + log(ℯ, (rij_1[i,j] / rij[i,j]))*sinij[i,j]
            )
            )
        end
        Vti[i] = Vti[i] + Vinf*(Panel_data.cos_theta[i]*cos(α) + Panel_data.sin_theta[i]*sin(α))                    #cos(Panel_data.theta[i] - α)
    end
    #println(Vti)

    return Vti
end
"""
    thetaij(
    panel_data
    )

Solves for the sin and cos of the difference between two angles for each panel relative to the other.

# Arguments:
- `panel_data::Panels` : struct with panel data

# Returns:
- `sinij::Array` : sin(theta_i - theta_j)
- `cosij::Array` : cos(theta_i - theta_j)
"""
function thetaij(
    panel_data
)
    n = length(panel_data.Panel_ID)
    sinij = Array{Float64, 2}(undef, n, n)
    cosij = Array{Float64, 2}(undef, n, n)
    for i = 1:n
        for j = 1:n
            sinij[i,j] = panel_data.sin_theta[i]*panel_data.cos_theta[j] - panel_data.cos_theta[i]*panel_data.sin_theta[j] #compute sin(theta_i - theta_j)
            cosij[i,j] = panel_data.cos_theta[i]*panel_data.cos_theta[j] + panel_data.sin_theta[i]*panel_data.sin_theta[j] #compute cos(theta_i - theta_j)
        end
    end
    return sinij, cosij
end

"""
    Coefficient_force_computer(
    Vti,
    Vinf,
    chord_length,
    panel_data
    )

Solves system of equations for the Hess-Smith Panel Method - see Computational Aerodynamics by Andrew Ning equation 2.234

# Arguments:
- `Vti::Vector` : Tangent velocity at each panel - see Tangent_velocity_computer
- `Vinf::Float` : Free stream velocity
- `chord_length::Float` : length of airfoil chord
- `panel_data::Panels` : struct with panel data

# Returns:
- `Cd::Vector` : 2d Drag coefficient
- `Cl::Vector` : 2d Lift coefficient
"""
function Coefficient_force_computer(
    Vti,
    Vinf,
    chord_length,
    panel_data
)
    #initialize vectors and variables
    Cpi = Vector{Float64}(undef, length(Vti))
    Pressure_i = [Vector{Float64}(undef, 2) for a in 1:length(Vti)]
    sumx = 0.0
    sumz = 0.0


    #compute the pressure and force per unit length of each panel
    for i = 1:length(Vti)
        Cpi[i] = 1 - (Vti[i] / Vinf)^2 #see equation 2.238 in Computational Aerodyanmics by Andrew Ning
        #println(Cpi[i])
        Pressure_i[i][1] = (Cpi[i])*sin(panel_data.theta[i])*panel_data.Panel_length[i]
        Pressure_i[i][2] = -(Cpi[i])*cos(panel_data.theta[i])*panel_data.Panel_length[i]
        sumx = sumx + Pressure_i[i][1]
        sumz = sumz + Pressure_i[i][2]
    end
    
    cd = sumx / chord_length
    cl = sumz / chord_length
    return cd, cl, Cpi
end

"""
    Hess_Smith_Panel(
    panel_data,
    Vinf,
    α,
    chord_length
    )

Solves for the 2d coefficient of lift and drag using the Hess-Smith Panel method

# Arguments:
- `panel_data::Panels` : struct with panel data
- `Vinf::Float` : Free stream velocity
- `α::Vector` : Angle of Attack (degrees)
- `chord_length::Float` : length of airfoil chord

# Returns:
- `Cd::Vector` : 2d Drag coefficient
- `Cl::Vector` : 2d Lift coefficient
"""
function Hess_Smith_Panel(
    panel_data,
    Vinf,
    α,
    chord_length
)
    Betaij = Beta_computer(panel_data)
    sinij, cosij = thetaij(panel_data)
    Aij, rij, rij_1 = Aij_computer(panel_data, Betaij, sinij, cosij)
    B = B_computer(panel_data, α, Vinf)
    solution = solve_system(Aij, B)
    Vti = Tangent_velocity_computer(panel_data, rij, rij_1, Betaij, sinij, cosij, solution, Vinf, α)
    cd, cl, Cpi = Coefficient_force_computer(Vti, Vinf, chord_length, panel_data)

    return cd, cl, Cpi
end

#Creates NACA coordinates using Airfoil AirfoilTools
Test1 = NACA4(2.0, 4.0, 12.0, false)
#x,z = naca4(Test1)

#call panel setup function
#=
test_panels = panel_setup(x,z, graph = true, graph_filename = "BasicProject\\TestGraph.png")
sinij, cosij = thetaij(test_panels)
cd, cl = Hess_Smith_Panel(test_panels, 1.0, 8.0, 1.0)
println(cd)
println(cl)
=#


# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = FLOWFoil.AirfoilTools.joukowsky(center, radius)

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = FLOWFoil.AirfoilTools.joukowsky_flow(
    center, radius, Vinf, alpha
)

# - Your Stuff - #
panel_joukowsky = panel_setup(x, y)
cd, cl, cp2 = Hess_Smith_Panel(panel_joukowsky, Vinf, alpha, radius)
push!(cp2, cp2[360])
# - Plot Stuff - #
pl = plot(; xlabel="x", ylabel="cp", yflip=true)
plot!(
    pl,
    x,
    surface_pressure_coefficient;
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
)
 plot!(x, cp2, label="Hess-Smith")
 savefig(pl, "BasicProject\\Hess_Smith_vs_Analytic_Solution.png")
println(" ")