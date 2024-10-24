using DelimitedFiles, ForwardDiff, Plots

import FLOWFoil.AirfoilTools

using .AirfoilTools
include("VortexPanel.jl") #include Vortex Panel (Hess Smith Panel method) for input data

##########BoundaryLayerFlow Overview##########
#=
The purpose of this program is to predict the boundary layer flow over an airfoil.
This is to be paired with VortexPannel.jl which will get the inviscid flow velocity

This will predict flow in two regions:
-Laminar: Use Thwait's method for incompressible laminar region.
-Turbulent: Use Head's method

Once θ for each control point (turbulent or laminar) then the normal flow equation will be re evaulated using
equation 3.37 in computational Aerodynamics

The panel method will then be re evaulated using the new normal flow boundary condition
Lift and drag values will then be predicted from these new panels.
=#

"""
    compute_laminar_delta(
    panel_data,
    Ve;
    ν = 0.0000148
    )

Solves for the laminar boundary layer thickness

# Arguments:
- `panel_data::NT` : Named tuple with panel data - note length measurements must be in meters
- `Ve::Vector` : Exit velocity or velocity at the edge of the boundary layer (m / s)

# Keyword Arguments:
- `ν::Float = 0.0000148` : Kinematic viscosity of fluid (cSt), automatically set to kinematic viscosity of air at 15 degrees celsius

# Returns:
- `δ_star::Vector` : Boundary layer thickness for each panel (assumed to be constant for the entire panel). Note that it will return 0's for panels in the turbulent flow region
- `dve_dx::Vector` : Exit velocity derivatives for each panel
- `θ::Vector` : Momentum thickness (meters) for each panel
- `H::Float` : H for the last panel in the laminar section, it is outputted to be used for the initial condition of the turbulent flow regime
"""
function compute_laminar_delta(
    panel_data,
    ve;
    ν = 0.0000148 #this is nu by the way not v, it is the kinematic viscosity of air (or fluid moving around body). It will change depending on the temperature and fluid
)

    #get reordered panel data
    panel_data_new = order_index_for_boundary_layer(panel_data, ve)   

    #initialize values and vectors
    n = length(panel_data_new.top_midpoints)
    dve_dx_top = Vector{Float64}(undef, n) .= 0.0
    dve_dx_bottom = similar(dve_dx_top) .= 0.0
    θ_top = similar(dve_dx_top, n) .= 0.0
    θ_bottom = similar(dve_dx_top, n) .= 0.0
    δ_star_top = similar(dve_dx_top, n) .= 0.0
    δ_star_bottom = similar(dve_dx_top, n) .= 0.0
    h = 0.0
    exit = false
    re_top = similar(θ_top) .= 0.0
    re_bottom = similar(θ_top) .= 0.0

    #redefine important parameters in named tuple
    important_parameters = (dve_dx_top = dve_dx_top, dve_dx_bottom = dve_dx_bottom, theta_top = θ_top, theta_bottom = θ_bottom,
     delta_star_top = δ_star_top, delta_star_bottom = δ_star_bottom, reynolds_top = re_top, reynolds_bottom = re_bottom
    )

    #Compute dve_dx for each panel
    #for now I will be taking simple derivatives
    for i = 1:n
        if i == 1
            #forward derivative
            important_parameters.reynolds_top[i] = abs((panel_data_new.exit_velocity_top[i]*
            panel_data_new.top_midpoints[i]) / (ν))
            h = panel_data_new.top_midpoints[i + 1] - panel_data_new.top_midpoints[i]
            important_parameters.dve_dx_top[i] = (panel_data_new.exit_velocity_top[i + 1] - panel_data_new.exit_velocity_top[i]) / h

            important_parameters.reynolds_bottom[i] = abs((panel_data_new.exit_velocity_bottom[i]*
            panel_data_new.bottom_midpoints[i]) / (ν))
            h = panel_data_new.bottom_midpoints[i + 1] - panel_data_new.bottom_midpoints[i]
            important_parameters.dve_dx_bottom[i] = (panel_data_new.exit_velocity_bottom[i + 1] - panel_data_new.exit_velocity_bottom[i]) / h
        elseif i < n && i > 1
            #central derivative
            important_parameters.reynolds_top[i] = abs((panel_data_new.exit_velocity_top[i]*
            panel_data_new.top_midpoints[i]) / (ν))
            h = panel_data_new.top_midpoints[i + 1] - panel_data_new.top_midpoints[i - 1]
            important_parameters.dve_dx_top[i] = (panel_data_new.exit_velocity_top[i + 1] - panel_data_new.exit_velocity_top[i - 1]) / h

            important_parameters.reynolds_bottom[i] = abs((panel_data_new.exit_velocity_bottom[i]*
            panel_data_new.bottom_midpoints[i]) / (ν))
            h = panel_data_new.bottom_midpoints[i + 1] - panel_data_new.bottom_midpoints[i - 1]
            important_parameters.dve_dx_bottom[i] = (panel_data_new.exit_velocity_bottom[i + 1] - panel_data_new.exit_velocity_bottom[i - 1]) / h

        else
            #backwards derivative
            important_parameters.reynolds_top[i] = abs((panel_data_new.exit_velocity_top[i]*
            panel_data_new.top_midpoints[i]) / (ν))
            h = panel_data_new.top_midpoints[i] - panel_data_new.top_midpoints[i - 1]
            important_parameters.dve_dx_top[i] = (panel_data_new.exit_velocity_top[i] - panel_data_new.exit_velocity_top[i - 1]) / h

            important_parameters.reynolds_bottom[i] = abs((panel_data_new.exit_velocity_bottom[i]*
            panel_data_new.bottom_midpoints[i]) / (ν))
            h = panel_data_new.bottom_midpoints[i] - panel_data_new.bottom_midpoints[i - 1]
            important_parameters.dve_dx_bottom[i] = (panel_data_new.exit_velocity_bottom[i] - panel_data_new.exit_velocity_bottom[i - 1]) / h
        end
    end
#=
    #compute initial condition using equation 3.100
    θ[1] = sqrt(0.075*ν / dve_dx[1])
    H = 0.0
    #solve equation 3.95 numerically for θ, likely using the classic Runga-Kutta method (fourth order Runga-Kutta)
    for i = 1:n - 1
        h = panel_data.panel_mid_points[i + 1][1] - panel_data.panel_mid_points[i][1]
        k1 = (0.225*ν) / (ve[i]*θ[i]) - (3*θ[i] / ve[i])*dve_dx[i]
        k2 = (0.225*ν) / (ve[i]*(θ[i] + 0.5*k1*h)) - (3*(θ[i] + 0.5*k1*h) / ve[i])*dve_dx[i]
        k3 = (0.225*ν) / (ve[i]*(θ[i] + 0.5*k2*h)) - (3*(θ[i] + 0.5*k2*h) / ve[i])*dve_dx[i]
        k4 = (0.225*ν) / (ve[i + 1]*(θ[i] + k3*h)) - (3*(θ[i] + k3*h) / ve[i + 1])*dve_dx[i + 1]
        θ[i + 1] = θ[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h

        #compute λ
        λ = ((θ[i]^2) / ν)*dve_dx[i] #equation 3.88
            #compute H using equations 3.96 and 3.97

        if exit == true
            λ = -0.2
            i = n - 1
        end

        if λ > 0. 1
            λ = 0.1
            exit = false
        elseif λ <= -0.1
            exit = true
        else
            exit = false
        end
        
        if exit == false
            if λ >= 0
                H = 2.61 - 3.75*λ-5.24*λ^2
            else
                H = 2.088 + (0.0731) / (0.14 + λ)
            end
            
            if 2.1 < H < 2.8 && abs(log10(re[i])) > -40.557 + 64.8066*H - 26.7538*H^2 + 3.3819*H^3 #check if transition occurs - equation 3.129
                δ_star[i] = 0.0
                exit = true
            else
                #compute δ* from H and θ (page 98 has the equation that relates the three variables)
                δ_star[i] = H*θ[i]
                if i == n - 1
                    λ = ((θ[i + 1]^2) / ν)*dve_dx[i + 1]
                    if λ >= 0
                        H = 2.61 - 3.75*λ-5.24*λ^2
                    else
                        H = 2.088 + (0.0731) / (0.14 + λ)
                    end
                    #output δ* for each panel
                    δ_star[i + 1] = H*θ[i + 1]
                end
            end

        else
            for j = i:n #set the remainder of the boundary layer thickness terms to 0
                δ_star[j] = 0.0
            end
        end
    end
    =#
    #return δ_star, dve_dx, θ, H
end

"""
    compute_turbulent_delta(
    panel_data,
    ve,
    δ_star,
    dve_dx,
    θ,
    H;
    ν = 0.0000148
    )

Solves for turbulent boundary layer thickness

# Arguments:
- `panel_data::NT` : Named tuple with panel data - note length measurements must be in meters
- `ve::Vector` : Exit velocity or velocity at the edge of the boundary layer (m / s)
- `δ_star::Vector` : boundary layer thickness values for laminar section with 0.0 given for sections with flow transition or seperation
- `dve_dx::Vector` : derivative of the exit velocity with respect to change in horizontal length (x) for each panel
- `θ::Vector` : momentum thickness (meters) for each panel
- `H::Float` : final H for the last panel in the laminar section - see equation 3.96 and 3.97

# Keyword Arguments:
- `ν::Float = 0.0000148` : kinematic viscosity of fluid (cSt), automatically set to kinematic viscosity of air at 15 degrees celsius

# Returns:
- `δ_star::Vector` : boundary layer thickness for each panel (assumed to be constant for the entire panel).
"""
function compute_turbulent_delta(
    panel_data,
    ve,
    δ_star,
    dve_dx,
    θ,
    H;
    ν = 0.0000148
)
    #determine number of panels which need boundary layer thickness
    n = length(δ_star)
    flow_seperated = false
    for i = 1:n
        #if δ_star[i] is 0 then it means it's in the turbulent regime and will be solved for
        #note that for now I am using a fourth order Runga Kutta method for this system of two ODE's
        if δ_star[i] == 0.0 && flow_seperated == false

            #solve for H1 equation 3.117
            if H <= 1.6
                H1 = 0.8234*(H-1.1)^(-1.287) + 3.3
            else
                H1 = 1.5501*(H-0.6778)^(-3.064) + 3.3
            end
            
            #solve equations 3.119 and 3.120 simultaneously using the Runga Kutta 4th order method for systems of ODE's
            h = panel_data.panel_mid_points[i + 1][1] - panel_data.panel_mid_points[i][1]
            cf = (0.246*10^(-0.678*H))*(abs((ve[i]*θ[i])) / ν)^(-0.268) #equation 3.122
            k11 = (0.0306 / θ[i])*(abs(H1 -3))^(-0.6169) - dve_dx[i]*(H1 / ve[i]) - ((cf / 2) - dve_dx[i]*(θ[i] / ve[i])*(H + 2))*(H1 / θ[i])
            k12 = (cf / 2) - dve_dx[i]*(θ[i] / ve[i])*(H + 2)
            H1_new = H1 + k11*h / 2
            θ_new = θ[i] + k12* h / 2
            k21 = (0.0306 / θ_new)*(abs(H1_new - 3))^(-0.6169) - dve_dx[i]*(H1_new / ve[i]) - (H1_new / θ_new)*((cf / 2) - dve_dx[i]*(θ_new / ve[i])*(H + 2))
            k22 = (cf / 2) - dve_dx[i]*(θ_new / ve[i])*(H + 2)
            H1_new = H1 + k21*h / 2
            θ_new = θ[i] + k22*h / 2
            k31 = (0.0306 / θ_new)*(abs(H1_new - 3))^(-0.6169) - dve_dx[i]*(H1_new / ve[i]) - (H1_new / θ_new)*((cf / 2) - dve_dx[i]*(θ_new / ve[i])*(H + 2))
            k32 = (cf / 2) - dve_dx[i]*(θ_new / ve[i])*(H + 2)
            H1_new = H1 + k31*h
            θ_new = θ[i] + k32*h
            k41 = (0.0306 / θ_new)*(abs(H1_new - 3))^(-0.6169) - dve_dx[i + 1]*(H1_new / ve[i + 1]) - (H1_new / θ_new)*((cf / 2) - dve_dx[i + 1]*(θ_new / ve[i + 1])*(H + 2))
            k42 = (cf / 2) - dve_dx[i + 1]*(θ_new / ve[i + 1])*(H + 2)
            
            H1 = H1 + (1/6)*(k11 + 2*(k21 + k31) + k41)*h
            θ[i + 1] = θ[i] + (1/6)*(k12 + 2*(k22 + k32) + k42)*h

            #check if H1 is < 3.3
            if H1 < 3.3
                println("Flow has seperated!") #if it is then flag that the flow might be seperating
                flow_seperated = true
            else
                println("Flow is fine")
            end

            #Solve for H using H1
            if H1 >= 5.3 && flow_seperated == false
                H = 0.86*(H1 - 3.3)^(-0.777) + 1.1
            elseif H1 < 5.3 && flow_seperated == false
                H = 1.1538*(H1 - 3.3)^(-0.326) + 0.6778
            else
                H = 0.0
            end

            #compute new δ_star using H = (δ_star / θ)
            δ_star[i + 1] = H*θ[i + 1]
        end
    end
    return δ_star
end

"""
    order_index_for_boundary_layer(
    panel_data,
    vti
    )

The purpose of this function is to sort the panels into two named tuples which contain the top and bottom coordinates of the
airfoil starting at the leading edge going to the trailing edge

# Arguments:
- `panel_data::NT` : Named tuple with panel data - note length measurements must be in meters
- `vti::Vector` : vector with the tangent velocity value for each panel

# Returns:
- `reordered_panel_data::NT` : named tuple with ordered panel information for the top and bottom surface of the airfoil
"""
function order_index_for_boundary_layer(
    panel_data,
    vti
)
    n = length(panel_data.panel_mid_points)
    front_edge_index = 0
    exit = false
    
    #find x value closest to 0
    min_value = 0.0
    for i = 1:n
        if i == 1
            min_value = panel_data.panel_mid_points[i][1]
        elseif abs(panel_data.panel_mid_points[i][1]) < min_value
            min_value = panel_data.panel_mid_points[i][1]
        else
            min_value = min_value
        end
    end

    #check which index hits 0
    for i = 1:n
        if panel_data.panel_mid_points[i][1] == min_value && exit == false
            front_edge_index = i
            exit = true
        end
    end
    #panel_data = (panel_ID = ID, panel_start_points = startpoints, panel_end_points = endpoints, panel_mid_points = midpoints,
    #panel_length = l, theta = θ, sin_theta = sin_t, cos_theta = cos_t)


    #initialize new top and bottom vectors
    panel_top_midpoints_new = Vector{Float64}(undef, front_edge_index) .= 0.0
    panel_bottom_midpoints_new = similar(panel_top_midpoints_new) .= 0.0
    ve_top_new = similar(panel_top_midpoints_new) .= 0.0
    ve_bottom_new = similar(panel_top_midpoints_new) .= 0.0

    #reverse order of vectors
    for i = 1:front_edge_index
        panel_top_midpoints_new[i] = panel_data.panel_mid_points[front_edge_index - i + 1][1]
        panel_bottom_midpoints_new[i] = panel_data.panel_mid_points[front_edge_index + i][1]
        ve_top_new[i] = vti[front_edge_index - i + 1]
        ve_bottom_new[i] = vti[front_edge_index + i]
    end

    reordered_panel_data = (top_midpoints = panel_top_midpoints_new, bottom_midpoints = panel_bottom_midpoints_new, exit_velocity_top = ve_top_new,
    exit_velocity_bottom = ve_bottom_new)
return reordered_panel_data
end

function plot_delta(
    panel_data,
    δ_star,
    filename
)
    n = length(panel_data.panel_mid_points)
    x = Array{Float64, 1}(undef, n)
    for i = 1:n
        x[i] = panel_data.panel_mid_points[i][1]
    end

    fig1 = plot(x, δ_star)
    savefig(fig1, filename) 
end

####### Main Function Calls ######
Test1 = NACA4(2.0, 4.0, 12.0, true)
x,z = naca4(Test1)
test_panels = panel_setup(x,z)
cd, cl, Cpi, vti = Hess_Smith_Panel(test_panels, 1.0, 0.0, 0.125)
#=δ_star, dvti_dx, θ, H = =# compute_laminar_delta(test_panels, vti)
#δ_star = compute_turbulent_delta(test_panels, vti, δ_star, dvti_dx, θ, H)
#plot_delta(test_panels, δ_star, "HessSmithProject\\delta_test.png")
#test = order_index_for_boundary_layer(test_panels, vti)

################################
println("")