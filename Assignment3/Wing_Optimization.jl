using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, Optim, SNOW

#Plan of Action
#The induced drag from VLM will be the function I am optimizing.
#Define a new function that calls VLM. It will be written in the
#form that SNOW can use
#Write all my constraints using equation 5.1 in the optimization textbook
#Set my f as my induced drag
#set lift force, pitch angle, wingspan, and speed as constraints

function Wing_Optimization!(g, x)
    CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area = VLM(x[1], x[2], x[3])

    g[1] = x[3] - 8.0  #wingspan has to be 8.0 meters. Basically x[2] - 8 has to be 0.
    g[2] = 0.5*1.007*wing_area*CFz - 1.7 #minimum lift must be 1.7 newtons.
    #this is just the lift equation using coefficient of lift. The density of air is 1.007 m3/kg for an alititude of 2000 meters
    return CDiff
end



function fx!(g, x) #This is my test function for optimization
    #SNOW wants functions stated like equation 5.1 in the optimization textbook on page 154

    #This function uses the optimizer to find the square root of 2. Notice that my constraint has to be less than or equal to 0
    # So I put the equation into the constraint and it optimizes the function. Then I just have the function itself be equal to 
    #the x value used.
    y = x[1]

    g[1] = x[1]^2 - 2
    return y
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

#working example for optim
#=e
f1(x) = x[1]^2 - 2*x[1]
x0 = [-1.0]
results = optimize(x->f1(first(x)), x0) 
min = minimum(results)
println(min)
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
    alpha = 5*pi/180 #set the angle of attack to 5 degrees

    Panels_span = 30
    Panels_chord = 15
    Spacing_type_span = Cosine()
    Spacing_type_chord = Uniform()
    Rref = [0.0,0.0,0.0]
    Vinf = 1.0
    wing_area = sum(panel_area)*2 #this gives wing area for lift calculations that I will use in the optimization
    ref = Reference(wing_area, mean(chord_distribution), span, Rref, Vinf)
    fs = Freestream(Vinf, alpha, beta, [0.0;0.0;0.0]) #Define freestream Parameters 

    #create the surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, Panels_span, Panels_chord; fc = fc, spacing_s = Spacing_type_span, spacing_c = Spacing_type_chord)
    surfaces = [surface]

    #perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric = true)
    CF, CM = body_forces(system; frame=Wind()) #compute near body forces
    CDiff = far_field_drag(system) #compute farfield drag
    CFx, CFy, CFz = CF
    CMx, CMy, CMz = CM

    return CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area #, grid, surface these I added as outputs for troubleshooting
end

#defining variables
pitch_angle = 5*pi/180 #set pitch angle to 5 degrees
min_lift = 1.7 #set the minimum lift weight to 1.7 Newtons.
wing_span = 8.0 #set the wingspan to 8 meters
speed = 1.0 #set the speed to 1 m/s
density = 1.007 #This is in kg/m3

x0 = [2.0, 8.0]
ng = 3
#xopt, fopt, info = minimize(Wing_Optimization!, x0, ng)

#This is a working snow example using the fx function
#=
x0 = [-0.5; -0.5]
ng = 2
xopt, fopt, info = minimize(fx!, x0, ng)
=#
CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area = VLM([0.0 0.0], [3.0 3.0], 8)
println("")
println(CFz)
println(wing_area)
#println(surface)