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
    CSVHeader = ["alpha" "Clx" "Cly" "Clz" "Cmx" "Cmy" "Cmz" "CDiff"] # don't add commas or else it will mess things up.
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

#plots lift to drag ratio or wing efficiency for the test
function LiftToDrag(CSVFile)
    Array_plot = readdlm(CSVFile, ',')
    n = length(Array_plot[:, 1]) - 1 #it is the the length of the first column minus the first entry which is a header
    Header = Array_plot[1, :]
    Alpha = Array_plot[2:n + 1, 1]
    clx = Array{Float64, 1}(undef, n) #initialize variables
    clz = Array{Float64, 1}(undef, n)
    CDiff = Array{Float64, 1}(undef, n)
    LiftToDrag = Array{Float64, 1}(undef, n)
    clx = Array_plot[2:n, 2]
    clz = Array_plot[2:n, 4]
    CDiff = Array_plot[2:n, 8]
    #Lift to drag ratio calculations
    for i = 1:length(clx)
        LiftToDrag[i] = clz[i]/(abs(clx[i])+abs(CDiff[i]))
        if LiftToDrag[i] > 5000 #I added in this condition because at 0 both coefficients should be 0 and so it is trying to divide 0 over 0
            LiftToDrag[i] = 0
        end
    end
    #=
    plot1 = plot(Alpha, LiftToDrag, label=false, xlabel="Alpha (degrees)", ylabel="Lift to Drag Ratio")
    #debugging
    #=
    println(length(LiftToDrag))
    println(length(Alpha))
    =# 
    savefig(plot1, "Assignment2\\$(filename)WingEfficiency.png")
    return plot1
    =#
    return LiftToDrag, Alpha
end
plot()
function ARTestRepeater(ARRange, wingspan, filename) #basically this function inputs a range of aspect ratios + wingspan and will output a graph with all the plots on it
    for i = 1:length(ARRange)
        CSV, CSVHeader= ARTest(ARRange[i], wingspan)
        WriteFile(CSV, "Assignment2\\AspectRatio_$(ARRange[i]).csv", CSVHeader)
        LTD, Alpha = LiftToDrag("Assignment2\\AspectRatio_$(ARRange[i]).csv")
        plot1 = plot!(Alpha, LTD, label= "Aspect Ratio: $(ARRange[i])", xlabel="Alpha(degrees)", ylabel="Lift to Drag Ratio")
        savefig(plot1, "Assignment2\\$(filename).png")
    end

end

#Testing functions
#CF, CM, CDiff = VLM(40, 2, 10, 1*pi/180, 0 , [0.0; 0.0; 0.0;])
#b = ARTest(2, 4)
#ACSVTest, CSVTestHeader = ARTest(5, 4)
#WriteFile(CSVTest, "Assignment2\\CSVTest1.csv", CSVTestHeader)
#plotter("Assignment2\\CSVTest1.csv", "AR5")
#=
println("X Direction")
println(CF[1, :])
println("Y Direction")
println(CF[2, :])
println("Z Direction")
println(CF[3, :])
=#
println("") # I added this so it doesn't print those huge vectors

ARTestRepeater([0.5 1 2 5], 4.0 , "AR_Analysis2")