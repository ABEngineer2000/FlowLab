using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles

#This function tests aspect ratio of a wing
function ARTest(AspectRatio, wingspan) #function assuming square wing
    c = AspectRatio * wingspan #solves for chord
    alpha = Array{Float16, 1}(undef, 30)
    CF = Array{Float64, 2}(undef, 3, 30)
    CM = Array{Float64, 2}(undef, 3, 30)
    CDiff = Array{Float64, 1}(undef, 30)
    alpha[1] = -15*pi/180 #start angle of attack at -15 degrees
    for i = 1:length(alpha)
        CF[1,i], CF[2,i], CF[3, i], CM[1, i], CM[2, i], CM[3, i], CDiff[i] = VLM(wingspan*c, c, wingspan, alpha[i], 0, [0.0; 0.0; 0.0]) #solve for lift, drag, and moment
        if i < 30 #I added this so it doesn't try to increment alpha past the number of elements in the array
            alpha[i + 1] = alpha[i] + 1*pi/180 #increment each angle of attack by 1 radian
        end
    end
    CSV = Array{Float64, 2}(undef, length(CF[1,:]), 8) # 8 columns for each component of force, moment as well as CDiff and alpha
    CSV[:,1] = alpha*180/pi #first column is angle of attack, order will go x y z force, then x y z moment, then CDiff
    CSV[:,2] = CF[1,:]
    CSV[:,3] = CF[2, :]
    CSV[:,4] = CF[3, :]
    CSV[:,5] = CM[1, :]
    CSV[:,6] = CM[2, :]
    CSV[:,7] = CM[3, :]
    CSV[:,8] = CDiff 
    CSVHeader = ["alpha" "Clx" "Cly" "Clz" "Cmx" "Cmy" "Cmz" "CDiff"]

    return CSV, CSVHeader
end

function VLM(Sref, cref, bref, alpha, beta, omega)
    xle = [0.0, 0.0] #first number is the position of the leading edge closest to the fuselage in the chordwise direction ,2nd is the same thing but with the leading edge at the wingtip
    yle = [0.0, bref/2] #spanwise direction
    zle = [0.0, 0.0] #vertical direction, use this for dihedral
    chord = [bref/2, bref/2] #first number is chord at the fueslage, the next is the chord at the wingtip.
    theta = [0.0, 0.0] #this is twist (rotation about y-axis) at the fueslage and wingtip respectively.
    phi = [0.0, 0.0] #This is rotation about the x-axis.
    #fc = fill((xc) -> 0, 2) # camberline function for each section, I don't think I need this

    Panels_span = 30
    Panels_chord = 10
    Spacing_type_span = Uniform()
    Spacing_type_chord = Uniform()
    Rref = [0,0,0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, Rref, Vinf)
    fs = Freestream(Vinf, alpha, beta, omega) #Define freestream Parameters

    #create the surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, Panels_span, Panels_chord, spacing_s = Spacing_type_span, spacing_c = Spacing_type_chord)
    surfaces = [surface]

    #perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric = true)
    CF, CM = body_forces(system) #compute near body forces
    CDiff = far_field_drag(system) #compute farfield drag
    CFx, CFy, CFz = CF
    CMx, CMy, CMz = CM

    return CFx, CFy, CFz, CMx, CMy, CMz, CDiff
end

function WriteFile(CSVArray, filename, CSVHeader) #write files to CSV
    open(filename, "w") do io
        writedlm(io, CSVHeader, ',')
        writedlm(io, CSVArray, ',')
    end
end

#need to finish this function
function plotter(CSVFile, filename, CSVHeader)
    Array_plot = readdlm(CSVFile, ',')
    n = size(Array_plot, 1) - 1
end

#Testing functions
#CF, CM, CDiff = VLM(40, 2, 10, 1*pi/180, 0 , [0.0; 0.0; 0.0;])
#b = ARTest(2, 4)
CSVTest, CSVTestHeader = ARTest(5, 4)
WriteFile(CSVTest, "Assignment2\\CSVTest1.csv", CSVTestHeader)

#=
println("X Direction")
println(CF[1, :])
println("Y Direction")
println(CF[2, :])
println("Z Direction")
println(CF[3, :])
=#
println("") # I added this so it doesn't print those huge vectors