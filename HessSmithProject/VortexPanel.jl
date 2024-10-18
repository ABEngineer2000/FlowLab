using Plots, LinearAlgebra, DelimitedFiles

import FLOWFoil.AirfoilTools

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
#=
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
=#
"""
    panel_setup(
    x,
    z;
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
- `panel_data::NT` : Named tuple with each panel ID, left, right, and midpoint coordinates. Also includes length coordinates and angle coordinates relative to the x axis.
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
    n = length(x)
    ID = Vector{Int64}(undef, n - 1)
    startpoints = [Vector{Float64}(undef, 2) for a in 1:n - 1]
    endpoints = similar(startpoints)
    midpoints = similar(startpoints)
    l = Vector{Float64}(undef, n - 1)
    θ = similar(l)
    sin_t = similar(l)
    cos_t = similar(l)

    #add data to the vectors that will create the Panels struct
    for i = 1:n - 1
        ID[i] = i
        startpoints[i] = [x[i], z[i]] # = vcat(startpoints, [[x[i],z[i]]])
        endpoints[i] = [x[i+1], z[i + 1]]# = vcat(endpoints, [[x[i + 1], z[i + 1]]])
        midpoints[i] = [(x[i] + x[i + 1])/2, (z[i] + z[i + 1])/2]#vcat(midpoints, [[(x[i] + x[i + 1])/2 , (z[i] + z[i + 1])/2]])
        l[i] = sqrt((x[i + 1] - x[i])^2 + (z[i + 1] - z[i])^2)
        θ[i] = atan(z[i + 1] - z[i], x[i + 1] - x[i])#push!(θ, atan(z[i + 1] - z[i], x[i + 1] - x[i]))
        sin_t[i] = (z[i + 1] - z[i]) / l[i]#push!(sin_t, (z[i + 1] - z[i]) / l[i + 1])
        cos_t[i] = (x[i + 1] - x[i]) / l[i]#push!(cos_t, (x[i + 1] - x[i]) / l[i + 1])
    end

    #Get rid of the 0's that I put in to initialize the vectors I'm using

    #create the struct
    panel_data = (panel_ID = ID, panel_start_points = startpoints, panel_end_points = endpoints, panel_mid_points = midpoints,
     panel_length = l, theta = θ, sin_theta = sin_t, cos_theta = cos_t)

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
    panel_data
    )

Computes Betaij values for an airfoil

# Arguments:
- `panel_data::NT` : Named tuple with data for the given airfoil

# Returns:
- `betaij::Array` : Array with βij values, see Computational Aerodynamics by Andrew Ning equation 2.212
"""
function Beta_computer(
    panel_data
)
    betaij = Array{Float64, 2}(undef, length(panel_data.panel_ID), length(panel_data.panel_ID)) #create n x n Matrix

    #using equation 2.212 this for loop computes βij
    for i = 1:length(panel_data.panel_ID)
        for j = 1:length(panel_data.panel_ID)
            if i == j
                betaij[i,j] = π*1
            else
                betaij[i,j] = atan(
                        ((panel_data.panel_start_points[j][1] - panel_data.panel_mid_points[i][1])*(panel_data.panel_end_points[j][2] - panel_data.panel_mid_points[i][2])
                         - (panel_data.panel_start_points[j][2] - panel_data.panel_mid_points[i][2])*(panel_data.panel_end_points[j][1] - panel_data.panel_mid_points[i][1]))
                         ,
                            ((panel_data.panel_start_points[j][1] - panel_data.panel_mid_points[i][1])*(panel_data.panel_end_points[j][1] - panel_data.panel_mid_points[i][1])
                            + (panel_data.panel_start_points[j][2] - panel_data.panel_mid_points[i][2])*(panel_data.panel_end_points[j][2] - panel_data.panel_mid_points[i][2]))                       
                    )
            end
        end
    end
    return betaij
end

"""
    Aij_computer(
    panel_data,
    betaij
    )

Computes Aij values from an airfoil

# Arguments:
- `panel_data::NT` : Named tuple with information for the given airfoil
- `betaij::Array` : Array with βij values, see Computational Aerodynamics by Andrew Ning equation 2.212
- `sinij::Array` : sinij Array - see thetaij
- `cosij::Array` : cosij Array - see thetaij

# Returns:
- `Aij::Matrix` : Array with Aij values, see Computational Aerodynamics by Andrew Ning 2.223, 2.233, and 2.234
"""
function Aij_computer(
    panel_data,
    betaij,
    sinij,
    cosij
    )
    n = length(panel_data.panel_ID)
    #initialize Aij Matrix to be an n + 1 x n + 1 matrix
    Aij = Matrix{Float64}(I, n + 1, n + 1)
    rij_1 = Matrix{Float64}(I, n, n)
    rij = Matrix{Float64}(I, n, n)
    #Computes n x n entries
    for i = 1:n
        for j = 1:n
            rij_1[i,j] = sqrt((panel_data.panel_end_points[j][1] - panel_data.panel_mid_points[i][1])^2 + (panel_data.panel_end_points[j][2] - panel_data.panel_mid_points[i][2])^2) #computes rij + 1
            rij[i,j] = sqrt((panel_data.panel_start_points[j][1] - panel_data.panel_mid_points[i][1])^2 + (panel_data.panel_start_points[j][2] - panel_data.panel_mid_points[i][2])^2) #computes rij
            Aij[i, j] = (
                log(ℯ, rij_1[i,j] / rij[i,j])*sinij[i,j] + betaij[i,j]*cosij[i,j]
                ) #computes Aij
            
        end
    end

    #computes [1 through n rows, n + 1 column] for Aij see equation 2.223
    for i = 1:n
        for j = 1:n
            Aij[i, n + 1] = Aij[i, n + 1] + log(ℯ, (rij_1[i,j] / rij[i,j]))*cosij[i,j] - betaij[i,j]*sinij[i,j]
        end
    end

    #computes [n + 1 row, 1 through n columns]
    for j = 1:n
        Aij[n + 1, j] = (
            betaij[1,j]*sinij[1,j] - log(ℯ, (rij_1[1,j] / rij[1,j]))*cosij[1,j]
            + betaij[n,j]*sinij[n,j] - log(ℯ, (rij_1[n,j] / rij[n,j]))*cosij[n,j]
            )
    end

    #computes [n + 1 row, n + 1 column]
    for j = 1:n
        Aij[n + 1, n + 1] = (
            Aij[n + 1, n + 1]
            + betaij[1,j]*cosij[1,j]
            + log(ℯ, (rij_1[1,j] / rij[1,j]))*sinij[1,j]
            ) + (
                betaij[n,j]*cosij[n,j]
                +  log(ℯ, (rij_1[n,j] / rij[n,j]))*sinij[n,j]
            )
    end

    return Aij, rij, rij_1
end


"""
    B_computer(
    panel_data,
    α,
    vinf
    )

Computes b values from an airfoil see Computational Aerodynamics by Andrew Ning equations 2.233 and 2.232

# Arguments:
- `panel_data::NT` : Named tuple with information for the given airfoil
- `α::Float` : User specified angle of attack in degrees
- `vinf::Float` : User specified free streamm velocity relative to the airfoil chord length

# Returns:
- `b::Matrix` : Matrix with B values
"""
function B_computer(
    panel_data,
    α,
    vinf,
    ) 
    α = α*π / 180 #convert α to radians
    n = length(panel_data.panel_ID)
    b = Matrix{Float64}(undef, n + 1, 1)
    
    #compute 1 through n rows
    for i = 1:n
        b[i,1] = 2*π*vinf*(panel_data.sin_theta[i]*cos(α) - panel_data.cos_theta[i]*sin(α))      #sin(Panel_data.theta[i] - α)
    end

    #compute n +1 row
    b[n + 1,1] = -2*π*Vinf*(panel_data.cos_theta[1]*cos(α) + panel_data.sin_theta[1]*sin(α) + panel_data.cos_theta[n]*cos(α) + panel_data.sin_theta[n]*sin(α))            #(cos(Panel_data.theta[1] - α) + cos(Panel_data.theta[n]  - α))
    
    return b
end

"""
    solve_system(
    Aij,
    b
    )

Solves system of equations for the Hess-Smith Panel Method - see Computational Aerodynamics by Andrew Ning equation 2.234

# Arguments:
- `Aij::Matrix` : Aij Matrix - see Aij_computer
- `b::Matrix` : B values - see B_computer

# Returns:
- `solution::Matrix` : n + 1 x 1 Matrix where 1 through n rows are the source strengths and the n + 1 row is the vortex strength
"""
function solve_system(
    Aij,
    b
)
    n = length(b)
    solution = Matrix{Float64}(undef, n, 1)
    solution = inv(Aij)*b
    return solution
end

"""
    Tangent_Velocity_computer(
    panel_data,
    Aij,
    βij,
    q_λ_vector,
    vinf,
    α
    )

Computes the tangent velocity at the control point of each of the panels

# Arguments:
- `panel_data::NT` : Named tuple with data on the airfoil and its panels
- `rij::Matrix` : rij - see Aij_computer
- `rij_1::Matrix` : rij_1 - see Aij_computer
- `βij::Matrix` : βij Matrix - see Beta_computer
- `sinij::Array` : sinij Array - see thetaij
- `cosij::Array` : cosij Array - see thetaij
- `q_λ_vector::Matrix` : Solution vector - see solve_system
- `vinf::Float` : Free stream velocity relative to the airfoil chord length
- `α::Float` : Angle of attack in degrees

# Returns:
- `vti::Vector` : n length Vector which holds the value of the tangent velocity at the control point for each panel
"""
function Tangent_velocity_computer(
    panel_data,
    rij,
    rij_1,
    βij,
    sinij,
    cosij,
    q_λ_vector,
    vinf,
    α
    )
    α = α*π/180 #convert α to radians
    n = length(panel_data.panel_ID)
    vti = Vector{Float64}(undef, n)

    #compute's tangent velocity at each of the panels see equation 2.237 in Computational Aerodynamics by Andrew Ning
    for i = 1:n
        for j = 1:n
            vti[i] = (vti[i] +
            (1 / (2*π))*q_λ_vector[j]*(
                βij[i,j]*sinij[i,j] - log(ℯ, (rij_1[i,j] / rij[i,j]))*cosij[i,j]
            )
            + (q_λ_vector[n + 1] / (2*π))*(
                βij[i,j]*cosij[i,j] + log(ℯ, (rij_1[i,j] / rij[i,j]))*sinij[i,j]
            )
            )
        end
        vti[i] = vti[i] + vinf*(panel_data.cos_theta[i]*cos(α) + panel_data.sin_theta[i]*sin(α))                    #cos(Panel_data.theta[i] - α)
    end
    #println(Vti)

    return vti
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
    n = length(panel_data.panel_ID)
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
    vti,
    vinf,
    chord_length,
    panel_data
    )

Solves system of equations for the Hess-Smith Panel Method - see Computational Aerodynamics by Andrew Ning equation 2.234

# Arguments:
- `vti::Vector` : Tangent velocity at each panel - see Tangent_velocity_computer
- `vinf::Float` : Free stream velocity
- `chord_length::Float` : Length of airfoil chord
- `panel_data::NT` : Named tuple with panel data

# Returns:
- `cd::Vector` : 2d Drag coefficient
- `cl::Vector` : 2d Lift coefficient
"""
function Coefficient_force_computer(
    vti,
    vinf,
    chord_length,
    panel_data
)
    #initialize vectors and variables
    cpi = Vector{Float64}(undef, length(vti))
    pressure_i = [Vector{Float64}(undef, 2) for a in 1:length(vti)]
    sumx = 0.0
    sumz = 0.0

    n = length(vti)
    #compute the pressure and force per unit length of each panel
    for i = 1:n
        cpi[i] = 1 - (vti[i] / vinf)^2 #see equation 2.238 in Computational Aerodyanmics by Andrew Ning
        #println(cpi[i])
        pressure_i[i][1] = (cpi[i])*sin(panel_data.theta[i])*panel_data.panel_length[i]
        pressure_i[i][2] = -(cpi[i])*cos(panel_data.theta[i])*panel_data.panel_length[i]
        sumx = sumx + pressure_i[i][1]
        sumz = sumz + pressure_i[i][2]
    end
    
    cd = sumx / chord_length
    cl = sumz / chord_length
    return cd, cl, cpi
end

"""
    Hess_Smith_Panel(
    panel_data,
    vinf,
    α,
    chord_length
    )

Solves for the 2d coefficient of lift and drag using the Hess-Smith Panel method

# Arguments:
- `panel_data::NT` : Named tuple with panel data
- `vinf::Float` : Free stream velocity
- `α::Vector` : Angle of Attack (degrees)
- `chord_length::Float` : length of airfoil chord

# Returns:
- `cd::Float` : 2d Drag coefficient
- `cl::Float` : 2d Lift coefficient
- `cpi::Vector` : Pressure coefficients for each panel
- `vti::Float` : Tangent velocity for each panel
"""
function Hess_Smith_Panel(
    panel_data,
    vinf,
    α,
    chord_length
)
    betaij = Beta_computer(panel_data)
    sinij, cosij = thetaij(panel_data)
    Aij, rij, rij_1 = Aij_computer(panel_data, betaij, sinij, cosij)
    b = B_computer(panel_data, α, vinf)
    solution = solve_system(Aij, b)
    vti = Tangent_velocity_computer(panel_data, rij, rij_1, betaij, sinij, cosij, solution, vinf, α)
    cd, cl, cpi = Coefficient_force_computer(vti, vinf, chord_length, panel_data)

    return cd, cl, cpi, vti
end
"""
    grid_convergence_study(
    thickness,
    max_camber,
    max_camber_position,
    initial_number_panels,
    max_number_panels,
    deltan;
    vinf = 1.0,
    α = 0.0,
    chord_length = 1.0,
    animate_convergence = false,
    graph_output = ""
    )

Performs a grid convergence study on a NACA airfoil.

# Arguments:
- `thickness::Float` : max thickness of the airfoil as a percentage of the chord length
- `max_camber::Float` : max camber as a percentage of the chord length
- `max_camber_position::Float` : max camber position as a percentage of the chord length
- `initial_number_panels::Int` : initial number of panels to be tested
- `max_number_panels::Int` : max number of panels before the function stops iterating
- `deltan::Int` : number of panels to add in between each iteration

# Keyword Arguments:
- `vinf::Float = 1.0` : value of the free stream velocity
- `α::Float = 0.0` : angle of attack in degrees
- `chord_length::Float = 1.0` : length of the chord
- `animate_convergence::Boolean = false` : if true then the function will output an animation to the desired file location
- `graph_output::String = ""` : desired file location for the output animation gif

# Returns:
- `cl::Array` : 2d Drag coefficient array for each amount of panels
- `cd::Array` : 2d Lift coefficient array for each amount of panels
- `n::Array` : number of panels for each lift and drag coefficient comptued
"""
function grid_convergence_study(
    thickness,
    max_camber,
    max_camber_position,
    initial_number_panels,
    max_number_panels,
    deltan;
    vinf = 1.0,
    α = 0.0,
    chord_length = 1.0,
    animate_convergence = false,
    graph_output = ""
)
    #initialize vectors
    maximum([0.0, 1.2])
    l_n = round(Int, (max_number_panels - initial_number_panels) / deltan) + 1 #value of the length of number of panels Array
    n = Array{Int64, 1}(undef, l_n) 
    cl = Array{Float64, 1}(undef, l_n)
    cd = similar(cl, l_n)
    #Compute cl and cd for each number of panels
    for i = 1:l_n
        if i == 1
            n[i] = initial_number_panels
        else
            n[i] = n[i - 1] + deltan
        panels = NACA4(max_camber, max_camber_position, thickness, false)
        x,z = naca4(panels, N = n[i])
        panel_init = panel_setup(x, z)
        cd[i], cl[i] = Hess_Smith_Panel(panel_init, vinf, α, chord_length)
        end
    end
    
    
    
    if animate_convergence == true
        anim = @animate for i = 1:l_n
         plot_cl = scatter!([n[i]], [cl[i]], markercolor = "blue", markersize = 4, legend = false, xlims = (0.0, n[l_n]), ylims = (0.0, 7.0), xlabel = "Number of Panels", ylabel = "Lift Coefficient", title = "NACA 1212 α = 4.0 degrees")
        end
        gif(anim, graph_output, fps = 15)
        plot() #reset the plot
    end
    
    return n, cl, cd
end

#n ,cl, cd = grid_convergence_study(12.0, 1.2, 4.0, 160, 5000, 100, animate_convergence = true, graph_output = "BasicProject\\5000_panel_animation.gif")

#println(n, cl, cd)

#Creates NACA coordinates using Airfoil AirfoilTools

############################ Main Function Calls#####################

Test1 = NACA4(2.0 , 4.0, 12.0, false)
x,z = naca4(Test1)

#call panel setup function

test_panels = panel_setup(x,z, graph = false, graph_filename = "BasicProject\\TestGraph.png")
sinij, cosij = thetaij(test_panels)
cd, cl, Cpi = Hess_Smith_Panel(test_panels, 1.0, 0.0, 1.0)
println(cl)
#=
push!(Cpi, Cpi[160])
plot2 = plot(x, Cpi)
savefig(plot2, "BasicProject\\symmetric.png")
println(cd)
println(cl)
=#
##########################################


#comparing with joukowsky_flow

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = AirfoilTools.joukowsky(center, radius)

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = AirfoilTools.joukowsky_flow(
    center, radius, alpha, Vinf
)

# - Your Stuff - #
panel_joukowsky = panel_setup(x, y)
cd, cl, cp2 = Hess_Smith_Panel(panel_joukowsky, Vinf, alpha, radius)
splice!(x, 1)
splice!(x, length(x))
splice!(x, length(x))
splice!(cp2, 1)
splice!(cp2, length(cp2))
splice!(surface_pressure_coefficient, 1)
splice!(surface_pressure_coefficient, length(surface_pressure_coefficient))
splice!(surface_pressure_coefficient, length(surface_pressure_coefficient))
println(length(x))
println(length(cp2))
cp2[342:350] .= 0.25
cp2[1:5] .= 0.25
# - Plot Stuff - #
pl = plot(; xlabel="x", ylabel="cp", yflip=true)
plot!(
    pl,
    x[1:342],
    surface_pressure_coefficient[1:342];
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
)
 plot!(x[1:342], cp2[1:342], label="Hess-Smith")
 savefig(pl, "HessSmithProject\\Hess_Smith_vs_Analytic_Solution2.png")
 
println(" ")