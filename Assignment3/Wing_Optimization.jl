using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, Optim, SNOW

#Plan of Action
#The induced drag from VLM will be the function I am optimizing.
#Define a new function that calls VLM. It will be written in the
#form that SNOW can use
#Write all my constraints using equation 5.1 in the optimization textbook
#Set my f as my induced drag
#set lift force, pitch angle, wingspan, and speed as constraints

function Wing_Optimization(g, x)
    CFx, CFy, CFz, CMx, CMy, CMz, CDiff = VLM(x[1], x[2])

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

#working example for optim
#=e
f1(x) = x[1]^2 - 2*x[1]
x0 = [-1.0]
results = optimize(x->f1(first(x)), x0) 
min = minimum(results)
println(min)
=#

function VLM(chord, span) #Performs a Vortex lattice analysis
    xle = [0.0, 0.0] #first number is the position of the leading edge closest to the fuselage in the chordwise direction ,2nd is the same thing but with the leading edge at the wingtip
    yle = [0.0, span/2] #spanwise direction
    zle = [0.0, 0.0] #vertical direction, use this for dihedral
    chordref = [chord, chord] #first number is chord at the fueslage, the next is the chord at the wingtip.
    theta = [0.0, 0.0] #this is twist (rotation about y-axis) at the fueslage and wingtip respectively.
    phi = [0.0, 0.0] #This is rotation about the x-axis.
    #fc = fill((xc) -> 0, 2) # camberline function for each section, I don't think I need this
    beta = 0.0
    alpha = 5*pi/180 #set the angle of attack to 5 degrees

    Panels_span = 30
    Panels_chord = 15
    Spacing_type_span = Cosine()
    Spacing_type_chord = Uniform()
    Rref = [0.0,0.0,0.0]
    Vinf = 1.0
    ref = Reference(chord*span, chord, span, Rref, Vinf)
    fs = Freestream(Vinf, alpha, beta, [0.0;0.0;0.0]) #Define freestream Parameters

    #create the surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, Panels_span, Panels_chord, spacing_s = Spacing_type_span, spacing_c = Spacing_type_chord)
    surfaces = [surface]

    #perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric = true)
    CF, CM = body_forces(system; frame=Wind()) #compute near body forces
    CDiff = far_field_drag(system) #compute farfield drag
    CFx, CFy, CFz = CF
    CMx, CMy, CMz = CM

    return CFx, CFy, CFz, CMx, CMy, CMz, CDiff
end

#defining variables
pitch_angle = 5*pi/180 #set pitch angle to 5 degrees
min_lift = 1.7 #set the minimum lift weight to 1.7 Newtons.
wing_span = 8.0 #set the wingspan to 8 meters
speed = 1.0 #set the speed to 1 m/s

#This is a working snow example using the fx function
#=
x0 = [-0.5; -0.5]
ng = 2
xopt, fopt, info = minimize(fx!, x0, ng)
=#

println("")