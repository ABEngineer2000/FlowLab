using DelimitedFiles, FiniteDifferences

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
    Ve;
    ν = 0.0000148 #this is nu by the way not v, it is the kinematic viscosity of air (or fluid moving around body). It will change depending on the temperature
)
    #initialize values and vectors
    n = length(Ve)
    dve_dx = Vector{Float64}(undef, n)
    θ = similar(dve_dx, n)
    δ_star = similar(dve_dx, n)

    #Compute dve_dx for each panel
    for i = 1:n
        if i == 1
            
        elseif i == n

        else

        end
    end
    #compute initial condition using equation 3.100
    θ_0 = sqrt(0.075*ν / dve_dx[1])
    
    #solve equation 3.95 numerically for θ, likely using the classic Runga-Kutta method

    #compute δ* from H and θ (page 98 has the equation that relates the three variables)

    #output δ* for each panel

    return δ_star
end


println("")