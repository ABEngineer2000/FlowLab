using DelimitedFiles, ForwardDiff

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

Solves for the 2d coefficient of lift and drag using the Hess-Smith Panel method

# Arguments:
- `panel_data::NT` : Named tuple with panel data - note length measurements must be in meters
- `Ve::Vector` : Exit velocity or velocity at the edge of the boundary layer (m / s)

# Keyword Arguments:
- `ν::Float = 0.0000148` : Kinematic viscosity of fluid (cSt), automatically set to kinematic viscosity of air at 15 degrees celsius

# Returns:
- `δ_star::Vector` : Boundary layer thickness for each panel (assumed to be constant for the entire panel). Note that it will return 0's for panels in the turbulent flow region
- `dve_dx::Vector` : Exit velocity derivatives for each panel
- `θ::Vector` : Boundary layer thickness (m) for each panel
- `H::Float` : H for the last panel in the laminar section, it is outputted to be used for the initial condition of the turbulent flow regime
"""
function compute_laminar_delta(
    panel_data,
    Ve;
    ν = 0.0000148 #this is nu by the way not v, it is the kinematic viscosity of air (or fluid moving around body). It will change depending on the temperature and fluid
)
    #initialize values and vectors
    n = length(Ve)
    dve_dx = Vector{Float64}(undef, n)
    θ = similar(dve_dx, n)
    δ_star = similar(dve_dx, n)
    h = 0.0
    exit = false

    #Compute dve_dx for each panel
    #for now I will be taking simple derivatives
    for i = 1:n
        if i == 1
            h = panel_data.panel_mid_points[i + 1][1] - panel_data.panel_mid_points[i][1]
            dve_dx[i] = (Ve[i + 1] - Ve[i]) / h #forward derivative
        elseif i < n && i > 1
            h = (panel_data.panel_mid_points[i + 1][1] - panel_data.panel_mid_points[i - 1][1])
           dve_dx[i] = 1*(Ve[i + 1] - Ve[i - 1]) / h  #central derivative
        else
            h = panel_data.panel_mid_points[i][1] - panel_data.panel_mid_points[i - 1][1]
           dve_dx[i] = (Ve[i] - Ve[i - 1]) / h #backwards derivative
        end
    end

    #compute initial condition using equation 3.100
    θ[1] = sqrt(0.075*ν / dve_dx[1])
    H = 0.0
    #solve equation 3.95 numerically for θ, likely using the classic Runga-Kutta method (fourth order Runga-Kutta)
    for i = 1:n - 1
        h = panel_data.panel_mid_points[i + 1][1] - panel_data.panel_mid_points[i][1]
        k1 = (0.225*ν) / (Ve[i]*θ[i]) - (3*θ[i] / Ve[i])*dve_dx[i]
        k2 = (0.225*ν) / (Ve[i]*(θ[i] + 0.5*k1*h)) - (3*(θ[i] + 0.5*k1*h) / Ve[i])*dve_dx[i]
        k3 = (0.225*ν) / (Ve[i]*(θ[i] + 0.5*k2*h)) - (3*(θ[i] + 0.5*k2*h) / Ve[i])*dve_dx[i]
        k4 = (0.225*ν) / (Ve[i]*(θ[i] + k3*h)) - (3*(θ[i] + k3*h) / Ve[i])*dve_dx[i]
        θ[i + 1] = θ[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h

        #compute λ
        λ = ((θ[i]^2) / ν)*dve_dx[i] #equation 3.88
        if λ <= -0.1
            exit = true
            for j = i:n
                δ_star[j] = 0.0
            end
            i = n - 1
        end
        if exit == false
            #compute H using equations 3.96 and 3.97
            if λ >= 0
                H = 2.61 - 3.75*λ-5.24*λ^2
            else
                H = 2.088 + (0.0731) / (0.14 + λ)
            end
            

            #compute δ* from H and θ (page 98 has the equation that relates the three variables)
                δ_star[i] = H*θ[i]
                if i == n - 1
                    λ = ((θ[i + 1]^2) / ν)*dve_dx[i + 1]
                    if λ <= -0.1
                        exit == true
                        for j = i:n
                            δ_star[j] = 0.0
                        end
                        i = n - 1
                    end
                    if exit == false
                        if λ >= 0
                            H = 2.61 - 3.75*λ-5.24*λ^2
                        else
                            H = 2.088 + (0.0731) / (0.14 + λ)
                        end
                        #output δ* for each panel
                        δ_star[i + 1] = H*θ[i + 1]
                        println(H)
                    end
                end
            end
    end
    return δ_star, dve_dx, θ, H
end


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
    for i = 1:n
        #if δ_star[i] is 0 then it means it's in the turbulent regime and will be solved for
        if δ_star[i] == 0.0
            #solve equation 3.119 and 3.120 for H1 and θ coupled
            
            #solve for H1 equation 3.117
            if H <= 1.6
                H1 = 0.8234*(H-1.1)^(-1.287) + 3.3
            else
                H1 = 1.5501*(H-0.6778)^(-3.064) + 3.3
            end
            
            cf = (0.246*10^(-0.678*H))*((ve[i]*θ[i]) / ν)^(-0.268) #equation 3.122
            k11 = (0.0306 / θ[i])*(H1 -3)^(-0.6169) - dve_dx[i]*(H1 / ve[i]) - ((cf / 2) - dve_dx[i]*(θ[i] / ve[i])*(H + 2))*(H1 / θ[i])
            k12 = (cf / 2) - dve_dx[i]*(θ[i] / ve[i])*(H + 2)
                #check if H1 is < 3.3
                
                #if it is then flag that the flow might be seperating

            #Solve for H using H1

            #solve equation 3.120 for θ

            #compute new δ_star using H = (δ_star / θ)
        end
    end

end

####### Main Function Calls ######
Test1 = NACA4(2.0, 4.0, 12.0, false)
x,z = naca4(Test1)
test_panels = panel_setup(x,z)
cd, cl, Cpi, Vti = Hess_Smith_Panel(test_panels, 1.0, 0.0, 1.0)
δ_star, dvti_dx, θ, H = compute_laminar_delta(test_panels, Vti)
println(H)
#compute_turbulent_delta(test_panels, Vti, δ_star, dvti_dx, θ, H)


################################

println("")