using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, SNOW, StaticArrays

#=
The plan for this project is to be able to design an aircraft that is both statically and
trim stable using the vortex lattice method, optimizing, and airfoil polars.

The idea is that first the vortex lattice analysis of the finite wing will be performed.
That data will be used to find the induced angle of attack for each cross section of the wing.
Once the induced angle of attack is calculated, the corresponding polar will be found using
a lookup table. The lift, moment, and drag calculations will then be performed on the aircraft.
The optimizer will take this data and optimize to make the aircraft have trim stability. 
The constraints will be a certain lift, drag, and static stability value to maintain.
The wing will be made of a wing, horizontal, and vertical stabilizer.

Design variables will be the chord lengths for the stabilizers.
=#

#doc string example

"""
    function call ie VLM(leading_edge_distribution, chord_distribution, twist_distribution, span, α, camber_line_function, S)

Function definition

**Arguments**
- `leading_edge_distribution::type`: A vector of leading edge values for each chord section

**Returns**
- `urfaces::type`: the surfaces output
- `system::type`: the system output
"""

function VLM(leading_edge_distribution, chord_distribution, twist_distribution, span, α, camber_line_function, S) #Performs a Vortex lattice analysis
    #this function requires two distribution vecotrs which contain the 2d wing geometry.
    xle = Array{Float64, 1}(undef, length(chord_distribution))
    yle = Array{Float64, 1}(undef, length(chord_distribution))
    zle = Array{Float64, 1}(undef, length(chord_distribution))
    chord = Array{Float64, 1}(undef, length(chord_distribution))
    theta = Array{Float64, 1}(undef, length(chord_distribution))
    phi = Array{Float64, 1}(undef, length(chord_distribution))
    panel_area = Array{Float64, 1}(undef, length(chord_distribution))
    for i = 1:length(chord_distribution)
        xle[i] = leading_edge_distribution[i] #each leading edge is according to the input leading edge distribution vector 
        yle[i] = (i - 1)*span/(length(chord_distribution)*2) #spanwise placement of the panels. Each planel is placed according to its top left corner
        zle[i] = 0.0 #vertical direction, I will assume a flat wing so my z coordinate is 0.
        chord[i] = chord_distribution[i] #first number is chord at the fueslage, the next is the chord at the wingtip.
        theta[i] = twist_distribution[i] #this is twist (rotation about y-axis) at the fueslage and wingtip respectively.
        phi[i] = 0.0 #This is rotation about the x-axis.
        panel_area[i] = chord_distribution[i] * span/(length(chord_distribution)*2) #outputs the area for each panel
    end

    fc = fill( camber_line_function, length(chord_distribution)) # camberline function for each section, it creates a camber in the z direction. make the function in terms of xc
    beta = 0.0
    alpha = α

    Panels_span = 50
    Panels_chord = 50
    Spacing_type_span = Uniform()
    Spacing_type_chord = Uniform()
    Rref = [0.0,0.0,0.0]
    Vinf = 1.0
    wing_area = S #0.11305 #this gives wing area for lift calculations that I will use in the optimization
    ref = Reference(wing_area, mean(chord_distribution), span, Rref, Vinf)
    fs = Freestream(Vinf, alpha, beta, [0.0;0.0;0.0]) #Define freestream Parameters 

    #create the surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, Panels_span, Panels_chord; fc = fc, mirror = false, spacing_s = Spacing_type_span, spacing_c = Spacing_type_chord)
    
    #println(grid)
    #println(surface)
    grids = [grid]
    surfaces = [surface]
    #rcp = VortexLattice.controlpoint(surfaces)

    #perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric = true)
    dCF, dCM = stability_derivatives(system)
    dCFz = dCF[1][1]
    dCMy = dCM[1][2]
    CF, CM = body_forces(system; frame=Wind()) #compute near body forces
    CDiff = far_field_drag(system) #compute farfield drag
    CFx, CFy, CFz = CF
    CMx, CMy, CMz = CM
    control_points_span = rcp_finder(surfaces, system, Panels_span, Panels_chord)
    Gamma = (system.Γ)[:]
    #println(length(Gamma))
    V_induced, α_i = Induced_AOA(control_points_span, surfaces, Gamma, Vinf)

    #testing getting the rcp values
    #=
    rcp1 = (system.surfaces[1])[2, 2].rcp
    println(rcp1)
    =#
    #Panel_Properties = get_surface_properties(system)
    #write_vtk("FinalProject\\symmetric-planar-wing", surfaces, Panel_Properties; symmetric = true)

    return #=surfaces, system, α_i,=# CFz, CMy, dCFz, dCMy #, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area, Panel_Properties #, grid, surface these I added as outputs for troubleshooting
end

function rcp_finder(surfaces, system, panels_span, panels_chord)
    #this inputs the Vortexlattice.jl surfaces and system and outputs a 3d matrix with the coordinates to each of the control points. It also
    #needs to know the number of panels spanwise and chordwise.
    #output will look like [control_point1;;; controlpoint2;;; controlpoint3;;; .... controlpointn]
    #rcp1 = (system.surfaces[1])[2, 2].rcp
    #println(rcp1)
    control_points_span = Array{Float64, 3}(undef, 1, 3, panels_span)
    n = 1 #this is added to index the control points
    #i goes along the span and j goes along the chord
    #option to take control points at chorter chord using this code in finding contro_points_span; round(Integer, panels_chord*0.25)
    for i = 1:panels_span
            control_points_span[1, :, n] = (system.surfaces[1])[round(Integer, panels_chord*0.40), i].rcp
            n = n + 1
    end
    return control_points_span
end

function Induced_AOA(control_points_span, surfaces, Γ, Vinf)
    #this function finds the induced angle of attack for each panel based off of the control points and circulation (Γ values).
    V_induced = Vector{SVector{3, Float64}}(undef, length(control_points_span[1, 1, :]))
    α_i = Array{Float64, 1}(undef, length(control_points_span[1, 1, :]))
    for i = 1:length(V_induced)
        #println(control_points[1, :, i])
        V_induced[i] = VortexLattice.induced_velocity(control_points_span[1, :, i], surfaces[1], Γ; symmetric = true)
        #V_induced[i] = VortexLattice.induced_velocity(i, surfaces[1], Γ; symmetric = true)
        α_i[i] = atan((V_induced[i])[3]/(Vinf))
        #println(V_induced[i][3])
        #println(α_i[i] * 180/pi + 2.0)
    end
    return V_induced, α_i
end

#This function will search the tabulated data for the matching lift coefficient
function Tabulated_Data_Match(AirfoilCSV, α_i)
    #input is a CSV file all the airfoil data as well as α_i which is the angle of attack in degrees
    data_cells, header_cells = readdlm(AirfoilCSV, ',', Float64, '\n'; header=true)
    α = data_cells[:,1]
    converged = data_cells[:, 5]
    amount_converged = 0
    #counts the number of converged indeces.
    for i = 1:length(α)
        if  converged[i] > 0.5
            amount_converged = amount_converged + 1
        else
            amount_converged = amount_converged
        end
    end
    #creates a new list of cl, cd, cm values that is the length of converged indeces.
    α_new = Array{Float64, 1}(undef, amount_converged)
    cl_new = Array{Float64, 1}(undef, amount_converged)
    cd_new = Array{Float64, 1}(undef, amount_converged)
    cm_new = Array{Float64, 1}(undef, amount_converged)
    j = 1 #new index I need to use for the following for loop
    #delete's all the non-converged data
    for i = 1:length(α)
        if converged[i] > 0.5
            α_new[j] = data_cells[i, 1]
            cl_new[j] = data_cells[i, 2]
            cd_new[j] = data_cells[i, 3]
            cm_new[j] = data_cells[i, 4]
            j = j + 1 #increment the new index
        else
            #nothing happens
        end
        #println(cl_new)
    end
    cl = Array{Float64, 1}(undef, length(α_i))
    cd = Array{Float64, 1}(undef, length(α_i))
    cm = Array{Float64, 1}(undef, length(α_i))
    
    #finds the closest value to the specified angle of attack
    #this is an example line that would find the index of the closest value of α to 3.05
    #min = findmin(broadcast(abs, (α .- 3.05)))[2]
    for i = 1:length(α_i)
    cl[i] = cl_new[findmin(broadcast(abs, (α_new .- α_i[i])))[2]]
    cd[i] = cd_new[findmin(broadcast(abs, (α_new .- α_i[i])))[2]]
    cm[i] = cm_new[findmin(broadcast(abs, (α_new .- α_i[i])))[2]]
    end

    return cl, cd, cm
end

#this function performs the wing analysis combining both the vortex lattice method for a finite wing and the vortex panel method for a 2d airfoil.
function Improved_wing_analysis(leading_edge_distribution, chord_distribution, span, AirfoilCSV, α) 
#this function inputs the leading edge distribution, chord distribution, span, and a CSV containing corresponding lift and drag coefficients for the 2d airfoil
#α is the angle of attack in radians
#calculates the wing area per section - note that this is the projected area of the section
section_area = Array{Float64, 1}(undef, length(chord_distribution))
section_span = span / length(chord_distribution)
for i = 1:length(section_area)
    section_area[i] = chord_distribution[i]*section_span
end
#call VLM to get the induced angle of attack for each section
surfaces, system, α_i, CFz, CMy, dCFz, dCMy = VLM(leading_edge_distribution, chord_distribution, span, α)

#find the lift coefficient for each section based on the induced angle of attack
#read in airfoil data, the first column is the angle of attack, the header will give the corresponding values of the columns
#this function assumes that the columns go as follows: alpha, c_l, c_d, c_m, and converged

cl_section = Array{Float64, 1}(undef, length(chord_distribution))
cd_section = Array{Float64, 1}(undef, length(chord_distribution))
cm_section = Array{Float64, 1}(undef, length(chord_distribution))
α_new = α_i .+ α
α_new = α_new .* 180/pi
cl_section, cd_section, cm_Section = Tabulated_Data_Match(AirfoilCSV, α_new)

#Calculate the new lift/drag/moment coefficient by taking a weighted average of each lift/drag/moment multiplied by their respective wing section area.
cl_sum = 0.0
cd_sum = 0.0
cm_sum = 0.0
wing_area = 0.0
for i = 1:length(chord_distribution)
    cl_sum = cl_sum + cl_section[i]*chord_distribution[i]
    cd_sum = cd_sum + cd_section[i]*chord_distribution[i]
    cm_sum = cm_sum + cm_section[i]*chord_distribution[i]
    wing_area = wing_area + section_area[i]
end
Cl = (cl_sum / sum(chord_distribution))*mean(chord_distribution) / wing_area
Cd = (cd_sum / sum(chord_distribution))*mean(chord_distribution) / wing_area
Cm = (cm_sum / sum(chord_distribution))*mean(chord_distribution) / wing_area

#output the new lift/drag coefficient that are based off the 3d definition
    return Cl, Cd, Cm, wing_area, CFz, CMy, dCFz, dCMy
end

#this function gathers a bunch of Cl data based on the range of angle of attack values specified
function GatherData(leading_edge_distribution, chord_distrubtion, span, AirfoilCSV, ComparisonData, α_range, α_resolution, save_file, CSVName)
    #all α values must be in radians
    #α_range must be a 2 element array with the starting and ending α values
    #creates all the α values that will be tested
    delta_α = (α_range[2] - α_range[1]) / α_resolution
    α = Array{Float64, 1}(undef, round(Integer, delta_α) + 1)
    for i = 1:length(α)
        if i < 2
            α[i] = α_range[1]
        else
            α[i] = α[i - 1] + α_resolution
        end
    end
    #Creates a Cl/Cd/Cm array to store all the Cl values
    Cl = Array{Float64, 1}(undef, length(α))
    Cd = Array{Float64, 1}(undef, length(α))
    Cm = Array{Float64, 1}(undef, length(α))
    wing_area = 0.0
    #call Improved_wing_analysis to calculate the Cl values
    for i = 1:length(α)
        Cl[i], Cd[i], Cm[i], wing_area = Improved_wing_analysis(leading_edge_distribution, chord_distribution, span, AirfoilCSV, α[i])
        println(α[i] * 180/pi)
    end
    #plot results in desired save file
    Cl_plot = plot(α .* 180/pi, Cl, label=false, xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient")
    Cl_CSV = Array{Float64, 2}(undef, length(α), 2)
    Cl_CSV[:, 1] = α .* 180/pi
    Cl_CSV[:, 2] = Cl
    open("$(CSVName).csv", "w") do io
        writedlm(io, Cl_CSV, ',')
    end
    savefig(Cl_plot, save_file)
end
#I'm just writing this function to plot the experimental data I got with the studies that I found
function PlotComparison(CSV1, CSV2, CSVLabels, filename)
    #Inputs both CSV files and then plots the subsequent results on top of each other
    #this assumes α is in degrees
    #pulls information from CSV files, assuming column 1 is α and column 2 is Cl
    Array1 = readdlm(CSV1, ',')
    Array2 = readdlm(CSV2, ',')
    α1 = Array1[:, 1]
    Cl1 = Array1[:, 2]
    α2 = Array2[:, 1]
    Cl2 = Array2[:, 2]
    #resets plot
    plot()
    plot1 = plot(α1, Cl1, label = CSVLabels[1], xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient")
    plot1 = plot!(α2, Cl2, label = CSVLabels[2], xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient")
    #save file based on the specified filename
    savefig(plot1, filename)   
end

#I didn't want to download another julia package so I just created my own fucntion to find the mean
function mean(x)
    sum = 0.0
    for i = 1:length(x)
        sum = sum + x[i]
    end
    m = sum/length(x)
    return m
end

function pitch_stability_analysis(dCFz_wing, dCMy_wing, x_wing,CL_wing, CMy_wing, dCFz_tail, dCMy_tail, x_tail,CL_tail, CMy_tail, tail_efficiency, tail_volume_coefficient)
    #inputs the dCF, dCM, and the distance (x/c) or the distance from the center of gravity to the mean aerodynamic center (quarter chord)
    #inputs the Cmac or moment coefficient about the aerodynamic center as well as the coefficient of lift to compute trim stability
    #the equations for stability are on pages 91-92 in the flight vehicle design textbook
    #also note that the tail_volume_coefficient is divided by x_tail because I multiply it in later again
    dCmg = -x_wing*dCFz_wing + dCMy_wing -x_tail*dCFz_tail + dCMy_tail #compute the static stability
    Cmg = CMy_wing + CMy_tail*tail_efficiency*tail_volume_coefficient - (x_wing*CL_wing) - (x_tail*CL_tail)*tail_efficiency*tail_volume_coefficient
    return dCmg, Cmg
end

function wing_area_calculator(chord_distribution, span)
    #computes the wing area using trapezoid interpolation
    wing_area = 0.0
    for i = 1:length(chord_distribution)
        if i < length(chord_distribution)
            wing_area = wing_area + ((chord_distribution[i] + chord_distribution[i + 1])/2)*(span/(length(chord_distribution) - 1))
        else
            wing_area = wing_area
        end
    end
    return wing_area
end
#computes the centroid of a wing given a chord distirubtion and span, assuming that the chords are interpolated using trapezoids
function mean_aerodynamic_center_calculator(leading_edge_distribution, chord_distribution, span)
    sum = 0.0
    wing_area = wing_area_calculator(chord_distribution, span)
    for i = 1:length(chord_distribution)
        if i < length(chord_distribution)
            #sum = area of trapezoid multiplied by the mean aerodynamic center of the trapezoidal wing
            sum = ((chord_distribution[i] + chord_distribution[i + 1])/2)*(span/(length(chord_distribution) - 1)) * ((mean([chord_distribution[i], chord_distribution[i + 1]]))*0.25 + leading_edge_distribution[i])
        else
            sum = sum
        end
    end
    mean_aerodynamic_center = sum / wing_area
    return mean_aerodynamic_center
end

#computes the center of gravity for a wing assuming that the wing sections are trapezoidal
function CG_calculator(chord_distribution, span, leading_edge_distribution)
    sum = 0.0
    wing_area = wing_area_calculator(chord_distribution, span)
    for i = 1:length(chord_distribution)
        #the formula for a trapezoid depends on knowing which side is shorter hence why I have to add an extra elseif statement
        if i < length(chord_distribution) && chord_distribution[i + 1] >= chord_distribution[i]
            #sum = area of trapezoid multiplied by the centroid of the trapezoid
            sum = ((chord_distribution[i] + chord_distribution[i + 1])/2)*(span/(length(chord_distribution) - 1))*(((chord_distribution[i] + 2*chord_distribution[i + 1])/(3*(chord_distribution[i] + chord_distribution[i + 1])))*(span/(length(chord_distribution) - 1)) + leading_edge_distribution[i])
        elseif i < length(chord_distribution) &&chord_distribution[i] > chord_distribution[i + 1]
            sum = sum = ((chord_distribution[i] + chord_distribution[i + 1])/2)*(span/(length(chord_distribution) - 1))*(((chord_distribution[i + 1] + 2*chord_distribution[i])/(3*(chord_distribution[i] + chord_distribution[i + 1])))*(span/(length(chord_distribution) - 1)) + leading_edge_distribution[i])
        else
            sum = sum
        end
    end
    CG = sum / wing_area
    #println(CG)
end
#CG_calculator([2.0, 1.0], 5.0, [0.0, 0.0])

#perform a stability optimization for pitch stability
function stability_optim(wing_chord_initial, wingspan_initial, tail_chord_initial, tail_span_initial, tail_distance_initial, wing_distance_initial,
    wing_density_thickness, tail_density_thickness, aircraft_weight_CG, lift_constraint, leading_edge_constraint, twist_constraint, air_density, airspeed; twist_initial = zeros(1, length(wing_chord_initial)), leading_edge_initial = zeros(1, length(wing_chord_initial)), tail_leading_edge_initial = zeros(1, length(tail_chord_initial)), α = 0.0*pi/180, camber_line_function = xc -> 0.0)
    #variables that the optimizer will change: wing chord lengths, wingspan, tail chord lengths, tail_span, twist (for each chord length), tail_distance, wing_distance
    #variables that will stay constant once inputted into the function: wing_density, tail_density, aircraft weight, lift_constraint, leading_edge_constraint
    #Aircraft_weight_CG is a vector containing the aircraft weight and center of gravity without the wings
    #define points in x0 where each set of variables starts
    global leading_edge_coordinate = 1
    global wing_chord_coordinate = length(leading_edge_initial) + 1
    global wingspan_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial)
    global tail_chord_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + 1
    global tailspan_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial)
    global twist_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + 1
    global tail_distance_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial)
    global wing_distance_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial)
    global α_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial) + length(α)
    global tail_leading_edge_coordinate = length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial) + length(α) + 1
    global Wing_density_thickness = wing_density_thickness
    global Tail_density_thickness = tail_density_thickness
    global Aircraft_weight_CG = aircraft_weight_CG
    global Lift_constraint = lift_constraint
    global Leading_edge_constraint = leading_edge_constraint
    global Twist_constraint = twist_constraint
    global Air_density = air_density
    global Airspeed = airspeed


    #set initial conditions vector
    x0 = Array{Float64, 1}(undef, length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial) + length(α) + length(tail_leading_edge_initial))
    x0[1:length(leading_edge_initial)] = leading_edge_initial
    x0[length(leading_edge_initial) + 1:length(leading_edge_initial) + length(wing_chord_initial)] = wing_chord_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial)] = wingspan_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + 1: length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial)] = tail_chord_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial)] = tail_span_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + 1: length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial)] = twist_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial)] = tail_distance_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial)] = wing_distance_initial
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial) + length(α)] = α
    x0[length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial) + length(α) + 1: length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(tail_distance_initial) + length(wing_distance_initial) + length(α) + length(tail_leading_edge_initial)] = tail_leading_edge_initial
    #x0 = [leading_edge_initial, wing_chord_initial, wingspan_initial, tail_chord_initial, tail_span_initial, twist_initial, α, tail_leading_edge_initial, wing_distance_initial, tail_distance_initial]
    #set up constraint and variable boundaries


    ng = 4
    lx = Array{Float64, 1}(undef, length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(α) + length(tail_leading_edge_initial) + length(wing_distance_initial) + length(tail_distance_initial))
    ux = Array{Float64, 1}(undef, length(leading_edge_initial) + length(wing_chord_initial) + length(wingspan_initial) + length(tail_chord_initial) + length(tail_span_initial) + length(twist_initial) + length(α) + length(tail_leading_edge_initial) + length(wing_distance_initial) + length(tail_distance_initial))
    lx[:] .= -Inf
    ux[:] .= Inf
    lg = Array{Float64, 1}(undef, 4)
    ug = Array{Float64, 1}(undef, 4)
    lg[:] .= -Inf
    ug[:] .= 0.0

    #set options up for the optimizer
    ip_options = Dict(
    "max_iter" => 300,
    "tol" => 1e-6
    )
    solver = IPOPT(ip_options)
    options = Options(;solver)
    
    xopt, fopt, info = minimize(stability_optim2!, x0, ng, lx, ux, lg, ug, options)
    return xopt, fopt, info
    
end
#testing stability_optim
#stability_optim([1 2 3],[1 2 3 4 5],1,1,1,1,1,1,1,[1 2 3],[1 2 3], [4 5], [6 7 8 9 10]; twist_initial = [1.0, 2.0, 3.0])

function stability_optim2!(g, x)
    #Constraint order: Lift constraint, leading_edge_constraint, then static stability constraint,
    wing_CG = CG_calculator(x[wing_chord_coordinate:wingspan_coordinate - 1], x[wingspan_coordinate], x[leading_edge_coordinate:wing_chord_coordinate - 1])
    tail_CG = CG_calculator(x[tail_chord_coordinate:tailspan_coordinate - 1], x[tailspan_coordinate], x[tail_leading_edge_coordinate:length(x)])
    tail_efficiency = 0.90
    wing_area = wing_area_calculator(x[wing_chord_coordinate:wingspan_coordinate - 1], x[wingspan_coordinate])
    CFz_wing, CMy_wing, dCFz_wing, dCMy_wing = VLM(x[leading_edge_coordinate:wing_chord_coordinate - 1], x[wing_chord_coordinate:wingspan_coordinate - 1], x[twist_coordinate:tail_distance_coordinate - 1], x[wingspan_coordinate], x[α_coordinate], camber_line_function, wing_area)
    tail_area = wing_area_calculator(x[tail_chord_coordinate:tailspan_coordinate - 1], x[tailspan_coordinate])
    CFz_tail, CMy_tail, dCFz_tail, dCMy_tail = VLM(x[tail_leading_edge_coordinate:length(x)], x[tail_chord_coordinate:tailspan_coordinate - 1], zeros(1, length(x[leading_edge_coordinate:length(x)])), x[tailspan_coordinate], x[α_coordinate], xc -> 0.0, tail_area)
    tail_volume_coefficient = tail_area / (wing_area * mean(x[wing_chord_coordinate:wingspan_coordinate - 1]))
    CG = (wing_CG * wing_area * Wing_density_thickness[1]* Wing_density_thickness[2] + tail_CG*tail_area*Tail_density_thickness[1]*Tail_density_thickness[2] + Aircraft_weight_CG[1]*Aircraft_weight_CG[2])/(Aircraft_weight_CG[1] + Wing_density_thickness[1] + Tail_density_thickness[1])
    x_wing = mean_aerodynamic_center_calculator(x[leading_edge_coordinate:wing_chord_coordinate - 1], x[wing_chord_coordinate:wingspan_coordinate - 1], x[wingspan_coordinate]) + x[wing_distance_coordinate] - CG
    x_tail = mean_aerodynamic_center_calculator(x[tail_leading_edge_coordinate:length(x)], x[tail_chord_coordinate:tailspan_coordinate - 1], x[tailspan_coordinate]) + x[tail_distance_coordinate] - CG
    dCmg, Cmg = pitch_stability_analysis(dCFz_wing, dCMy_wing, x_wing, CFz_wing, CMy_wing, dCFz_tail, dCMy_tail, x_tail, CFz_tail, CMy_tail, tail_efficiency, tail_volume_coefficient)
    
    g[1] = dCmg #static stability constraint
    #lift constraint
    g[2] = Lift_constraint - (CFz_wing*wing_area + CFz_tail*tail_area)*0.5*Air_density*Airspeed^2
    #leading_edge_constraint\
    for i = 3:length(leading_edge_constraint) + 3
        g[i] = leading_edge_constraint - x[leading_edge_coordinate + i - 3] #Leading_edge_constraint - x[leading_edge_coordinate:wing_chord_coordinate - 1]
    end
    #twist_constraint
    for i = length(leading_edge_constraint) + 4:length(leading_edge_constraint) + length(twist_constraint)
        g[i] = Twist_constraint - x[twist_coordinate + i - length(leading_edge_constraint) - 4]
    end
    CMG = Cmg^2
    return CMG
end

leading_edge_distribution = Array{Float64, 1}(undef, 2)
chord_distribution = Array{Float64, 1}(undef, 2)
twist_distribution = Array{Float64, 1}(undef, 2)
leading_edge_distribution[:] .= 0.0
chord_distribution[:] .= 0.190 + 0.0475
twist_distribution[:] .= 0.0
span = 0.595
HS_distribution = Array{Float64, 1}(undef, 2)
HS_chord_distribution = Array{Float64, 1}(undef, 2)
HS_twist_distribution = Array{Float64, 1}(undef, 2)
HS_twist_distribution .= 0.0
HS_distribution[:] .= 0.0
HS_chord_distribution[:] .= 0.02
HS_span = 0.125
HS_location = 0.15
M = 0.06
p = 0.4
camber_line_function = xc -> begin xc < p ? (M/p^2)*(2*p*xc - xc^2) : (M/(1-p)^2)*(1- 2*p + 2*p*xc - xc^2) end  
camber_line_function_tail = xc -> 0.0
wing_area = wing_area_calculator(chord_distribution, span)
tail_area = wing_area_calculator(HS_chord_distribution, HS_span)
tail_volume_coefficient = tail_area / (wing_area*mean(chord_distribution))
wing_location = -0.8
lift_constraint = 20.0
leading_edge_constraint = 0.0
twist_constraint = 0.0
air_density = 0.03
airspeed = 30.0

stability_optim(chord_distribution, span, HS_chord_distribution, HS_span, HS_location, wing_location, [1.0, 1.0], [1.0, 1.0], [20.0, 3.0], lift_constraint, leading_edge_constraint, twist_constraint, air_density, airspeed, α = 2.0*pi/180)
#testing stability_optim2!
#stability_optim2!(1, [leading_edge_distribution, chord_distribution, span, HS_chord_distribution, HS_span, HS_location, wing_location, twist_distribution, 2.0*pi/180, leading_edge_distribution], camber_line_function)

#surfaces, system, α_i = VLM(leading_edge_distribution, chord_distribution, span)
#Cl, Cd, Cm, wing_area, CFz, CMy, dCFz, dCMy = Improved_wing_analysis(leading_edge_distribution, chord_distribution, span, "FinalProject\\Tabulated_Airfoil_Data\\NACA_6412.csv", 2.0*pi/180)
#=
CFz_wing, CMy_wing, dCFz_wing, dCMy_wing = VLM(leading_edge_distribution, chord_distribution, twist_distribution, span, 2.0*pi/180, camber_line_function, wing_area)
CFz_tail, CMy_tail, dCFz_tail, dCMy_tail = VLM(HS_distribution, HS_chord_distribution, HS_twist_distribution, HS_span, 2.0*pi/180, camber_line_function_tail, tail_area)
=#
#println(CFz_wing)
#println(Cl)
#println(CMy)
#GatherData(leading_edge_distribution, chord_distribution, span,"FinalProject\\Tabulated_Airfoil_Data\\NACA_6412.csv", ComparisonData, [-6.0*pi/180 15.0*pi/180], 1*pi/180, "FinalProject\\Accuracy_Comparison\\ComparitiveStudy_Comparisondata.png", "FinalProject\\Accuracy_Comparison\\ComparitiveStudy_Comparisondata")
#PlotComparison("FinalProject\\Accuracy_Comparison\\ComparitiveStudyFigure23.csv", "FinalProject\\Accuracy_Comparison\\ComparitiveStudy_Comparisondata.csv", ["Krishnan et al Data", "Improved Wing Analysis Data"], "FinalProject\\Accuracy_Comparison\\ComparisonStudy_vs_ImprovedWingAnalysis.png")
#this is to test my camber line function
#=
fc = fill((xc) -> begin xc < p ? (M/p^2)*(2*p*xc - xc^2) : (M/(1-p)^2)*(1- 2*p + 2*p*xc - xc^2) end, length(chord_distribution))
println(fc[1](1.0))
=#
#=
static_stability, trim_stability = pitch_stability_analysis(dCFz_wing, dCMy_wing, -0.8, CFz_wing, CMy_wing, dCFz_tail, dCMy_tail, 0.15, CFz_tail, CMy_tail, 0.90, tail_volume_coefficient)
println(static_stability)
println(trim_stability)
=#

lx = -Inf
ux = Inf
lg = -Inf
ug = 1.1
#this was created for some IPOPT tests I performed

#=
function fx1!(g, x)
    y = x[1] + x[2]
    g[1] = a*x[1]^2 - 2*x[2][1]
    return y
end
options = Options(solver=IPOPT())  # choosing IPOPT solver
a = 0.3
x0 = [2.1, 2.3]
xopt, fopt, info = minimize(fx1!, x0, 2)
=#

println("Done") #this is so I don't print anything I don't want.