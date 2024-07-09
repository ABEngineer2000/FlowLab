using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, SNOW, StaticArrays, FindClosest

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

function VLM(leading_edge_distribution, chord_distribution, span) #Performs a Vortex lattice analysis
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
        yle[i] = (i - 1)*span/length(chord_distribution) #spanwise placement of the panels. Each planel is placed according to its top left corner
        zle[i] = 0.0 #vertical direction, I will assume a flat wing so my z coordinate is 0.
        chord[i] = chord_distribution[i] + 0.0 #first number is chord at the fueslage, the next is the chord at the wingtip.
        theta[i] = 0.0 #this is twist (rotation about y-axis) at the fueslage and wingtip respectively.
        phi[i] = 0.0 #This is rotation about the x-axis.
        panel_area[i] = chord_distribution[i] * span/(length(chord_distribution) * 2) #outputs the area for each panel
    end
    fc = fill((xc) -> 0, length(chord_distribution)) # camberline function for each section, it creates a camber in the z direction. make the function in terms of xc
    beta = 0.0
    alpha = 2*pi/180 #set the angle of attack to 5 degrees

    Panels_span = length(xle)
    Panels_chord = 3
    Spacing_type_span = Cosine()
    Spacing_type_chord = Uniform()
    Rref = [0.0,0.0,0.0]
    Vinf = 1.0
    wing_area = sum(panel_area)*2 #this gives wing area for lift calculations that I will use in the optimization
    ref = Reference(wing_area, mean(chord_distribution), span, Rref, Vinf)
    fs = Freestream(Vinf, alpha, beta, [0.0;0.0;0.0]) #Define freestream Parameters 

    #create the surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, Panels_span, Panels_chord; fc = fc, spacing_s = Spacing_type_span, spacing_c = Spacing_type_chord)
    
    #println(grid)
    println()
    #println(surface)

    surfaces = [surface]
    #rcp = VortexLattice.controlpoint(surfaces)

    #perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric = true)
    CF, CM = body_forces(system; frame=Wind()) #compute near body forces
    CDiff = far_field_drag(system) #compute farfield drag
    CFx, CFy, CFz = CF
    CMx, CMy, CMz = CM
    control_points_span = rcp_finder(surfaces, system, Panels_span, Panels_chord)
    Gamma = (system.Γ)[:]
    V_induced, α_i = Induced_AOA(control_points_span, surfaces, Gamma, Vinf)

    #testing getting the rcp values
    #=
    rcp1 = (system.surfaces[1])[2, 2].rcp
    println(rcp1)
    =#
    #Panel_Properties = get_surface_properties(system)
    #write_vtk("FinalProject\\symmetric-planar-wing", surfaces, Panel_Properties; symmetric = true)

    return surfaces, system, α_i #CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area, Panel_Properties #, grid, surface these I added as outputs for troubleshooting
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
    for i = 1:panels_span
            control_points_span[1, :, n] = (system.surfaces[1])[1, i].rcp
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
        V_induced[i] = VortexLattice.induced_velocity(control_points_span[1, :, i], surfaces[1], Γ)
        α_i[i] = atand((V_induced[i])[3]/Vinf)
        #println(α_i[i])
    end
    return V_induced, α_i
end

#this function performs the wing analysis combining both the vortex lattice method for a finite wing and the vortex panel method for a 2d airfoil.
function Improved_wing_analysis(leading_edge_distribution, chord_distribution, span, AirfoilCSV) 
#this function inputs the leading edge distribution, chord distribution, span, and a CSV containing corresponding lift and drag coefficients for the 2d airfoil

#calculates the wing area per section - note that this is the projected area of the section
section_area = Array{Float64, 1}(undef, length(chord_distribution))
section_span = span / length(chord_distribution)
for i = 1:length(section_area)
    section_area[i] = chord_distribution[i]*section_span
end
#call VLM to get the induced angle of attack for each section
surfaces, system, α_i = VLM(leading_edge_distribution, chord_distribution, span)

#find the lift coefficient for each section based on the induced angle of attack
#read in airfoil data, the first column is the angle of attack, the header will give the corresponding values of the columns
#this function assumes that the columns go as follows: alpha, c_l, c_d, c_m, and converged
data_cells, header_cells = readdlm(AirfoilCSV, ',', Float64, '\n'; header=true)
α = data_cells[:,1]
cl_section = Array{Float64, 1}(undef, length(chord_distribution))
cd_section = Array{Float64, 1}(undef, length(chord_distribution))
cm_section = Array{Float64, 1}(undef, length(chord_distribution))

#this is an example line that would find the index of the closest value of α to 3.05
#min = findmin(broadcast(abs, (α .- 3.05)))[2]

for i = 1:length(chord_distribution)
    cl_section[i] = data_cells[findmin(broadcast(abs, (α .- α_i[i])))[2], 2]
    cd_section[i] = data_cells[findmin(broadcast(abs, (α .- α_i[i])))[2], 3]
    cm_section[i] = data_cells[findmin(broadcast(abs, (α .- α_i[i])))[2], 4]
end
#Calculate the new lift/drag coefficient by adding all the lift/drag coefficients together and dividing by the all the wing section areas added together
#output the new lift/drag coefficient

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

leading_edge_distribution = [0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0]
chord_distribution = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
span = 10.0
#surfaces, system, α_i = VLM(leading_edge_distribution, chord_distribution, span)
Improved_wing_analysis(leading_edge_distribution, chord_distribution, span, "FinalProject\\Tabulated_Airfoil_Data\\NACA_6412.csv")

println("Done") #this is so I don't print anything I don't want.