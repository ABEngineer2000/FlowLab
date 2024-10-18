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


function compute_laminar_delta(
    x,
    Ve;
    ν = 0.0000148 #this is nu by the way not v, it is the kinematic viscosity of air (or fluid moving around body). It will change depending on the temperature
)
    #initialize values and vectors
    n = length(Ve)
    dve_dx = Vector{Float64}(undef, n)
    θ = similar(dve_dx, n)
    δ_star = similar(dve_dx, n)

    #Compute dve_dx for each panel
    #for now I will be taking simple derivatives
    for i = 1:n
        if i == 1
            h = x.Panel_mid_points[i + 1][1] - x.Panel_mid_points[i][1]
            dve_dx[i] = (Ve[i + 1] - Ve[i]) / h #forward derivative
        elseif i == n
            h = (x.Panel_mid_points[i + 1] - x.Panel_mid_points[i - 1])
            dve_dx[i] = (Ve[i + 1] - Ve[i - 1]) / h
            #central derivative
        else
            h = x.Panel_mid_points[i] - x.Panel_mid_points[i - 1]
            dve_dx[i] = (Ve[i] - Ve[i]) / h
        end
    end

    #compute initial condition using equation 3.100
    θ[1] = sqrt(0.075*ν / dve_dx[1])
    
    #solve equation 3.95 numerically for θ, likely using the classic Runga-Kutta method
    for i = 1:n - 1
        h = x.Panel_mid_points[i + 1][1] - x.Panel_mid_points[i][1]
        k1 = (0.225*ν) / (Ve[i]*θ[i]) - (3*θ_0 / Ve[i])*dve_dx[i]
        k2 = (0.225*ν) / (Ve[i]*(θ[i] + 0.5*k1*h)) - (3*(θ[i] + 0.5*k1*h) / Ve[i])*dVe_dx[i]
        k3 = (0.225*ν) / (Ve[i]*(θ[i] + 0.5*k2*h)) - (3*(θ[i] + 0.5*k2*h) / Ve[i])*dVe_dx[i]
        k4 = (0.225*ν) / (Ve[i]*(θ[i] + k3*h)) - (3*(θ[i] + k3*h) / Ve[i])*dVe_dx[i]
        θ[i + 1] = θ[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h
    end

    #compute δ* from H and θ (page 98 has the equation that relates the three variables)

    #output δ* for each panel

    return δ_star
end



println("")