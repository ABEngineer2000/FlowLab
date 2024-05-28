using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, FiniteDiff

function Tail_Analysis(chord_length, wingspan, HS_chord, HS_span, VS_chord, VS_span, tail_distance, wing_distance, n, alpha_range, filename) #Assuming you have a rectangular wing and tail. HS 
    #is horizontal stabilizer whilst VS is vertical stabilizer. Tail distance is the distance from the center
    #of lift of the tail to the center of gravity whilst wing distance is the distance from the center
    #of lift of the wing to the center of gravity
    #n is the tail efficiency, see equation 5.23 in Flight Vehicle Design. Typical values range from 0.85-0.9
        #call VLM to get all the coefficients for wing
        #call VLM to get all the coefficients for HS
        #call VLM to get all the coefficients for VS
        #Compute numerical derivative's: dclz/dalpha, dcmx/dalpha, dcmy/dalpha, and if needed the roll stability: dclmz/dalpha.
        #output these derivatives
        #can also output an optional stability value for pitch. Bascically it will give dCcg/dalpha.
        #This uses equation 5.23 in the Flight Vehicle Design textbook
        #lastly it will output the tail volume coefficients for the horizontal and vertical stabilizer using Equation 5.33 in the Flight Vehicle Design Textbook
        wCFx, wCFy, wCFz, wCMx, wCMy, wCMz, wCDiff = VLM_AlphaRange(wingspan, chord_length, alpha_range) #w is for wing
        hCFx, hCFy, hCFz, hCMx, hCMy, hCMz, hCDiff = VLM_AlphaRange(HS_span, HS_chord, alpha_range) #h is for horizontal stabilizer
        vCFx, vCFy, vCFz, vCMx, vCMy, vCMz, vCDiff = VLM_AlphaRange(VS_span, VS_chord, alpha_range) # v is for vertical stabilizer

        #computes stability derivatives
        wing_dclz_dalpha = Derivative_calc(alpha_range, wCFz)
        wing_dcmy_dalpha = Derivative_calc(alpha_range, wCMy)

        HS_dclz_dalpha = Derivative_calc(alpha_range, hCFz)
        HS_dcmy_dalpha = Derivative_calc(alpha_range, hCMy)

        VS_dclz_dalpha = Derivative_calc(alpha_range, vCFz) #don't need either of these for pitch stability
        VS_dcmy_dalpha = Derivative_calc(alpha_range, vCMy)

        #compute tail volume coefficient. See page 94 Flight Vehicle Design
        Tail_Volume_Coefficient = (HS_span * HS_chord)/(wingspan * chord_length) #this is based on pitch stability

        #computes trim stability
        Cmg = Pitch_Trim_Stability(wCMy, wCFz, hCMy, hCFz, wing_distance, tail_distance, alpha_range, Tail_Volume_Coefficient,n, HS_chord)

        #compute stability derivatives. Note that this is assuming the cmy moments oppose the moments created by clz
        DCmg_Dalpha = Array{Float64, 1}(undef, length(wing_dclz_dalpha))
        for i = 1:length(wing_dclz_dalpha)
            DCmg_Dalpha[i] = wing_dcmy_dalpha[i] + HS_dcmy_dalpha[i]  - (wing_distance/chord_length)*wing_dclz_dalpha[i] - n*Tail_Volume_Coefficient*HS_dclz_dalpha[i]
        end
        #creates a CSV with all the stability data
        CSV = Array{Float64, 2}(undef, length(DCmg_Dalpha), 2)
        CSV[:, 1] = alpha_range
        CSV[:, 2] = DCmg_Dalpha
        CSV_Header = ["Alpha (Degrees)" "DCmg_Dalpha"]
        WriteFile(CSV, filename, CSV_Header)
        return Tail_Volume_Coefficient, DCmg_Dalpha, Cmg #returns tail volume coefficient and main stability derivative
end

function Tail_Volume_Coefficient_Repeater(chord_length, wingspan, HS_chord, HS_span, VS_chord, VS_span, tail_distance, wing_distance, n, alpha_range, filename, Tail_Volume_Scalars) #This function will call Tail_Analysis for ever tail volume coefficient you want to analyze
    #after that, it will plot each solution on a curve and save that curve to desired file destination
    #it will specifically change the tail volume coefficient by multiplying both horizontal tail and chord by all the tail volume scalars
    #it will calculate the stability derivatives about the center of gravity
    #note that alpha range must include all values of alpha you want to analyze
    filename1 = Array{String, 1}(undef, length(Tail_Volume_Scalars))
    HS_span1 = Array{Float64, 1}(undef, length(Tail_Volume_Scalars))
    HS_chord1 = Array{Float64, 1}(undef, length(Tail_Volume_Scalars))
    Tail_Volume_Coefficient = Array{Float64, 1}(undef, length(Tail_Volume_Scalars))
    DCmg_Dalpha = Array{Float64, 2}(undef, length(alpha_range), length(Tail_Volume_Scalars))
    Cmg = Array{Float64, 2}(undef, length(alpha_range), length(Tail_Volume_Scalars)) #This is for trim stability
    plot() #resets plot
    for i =1:length(Tail_Volume_Scalars)
        HS_span1[i] = HS_span * Tail_Volume_Scalars[i] #Creates new Horizontal Stabilizer coefficients based on the Tail Volume Scalars
        HS_chord1[i] = HS_chord * Tail_Volume_Scalars[i]
        Tail_Volume_Coefficient[i] = (HS_span1[i] * HS_chord1[i])/(wingspan * chord_length) #Computes new tail volume coefficient
        filename1[i] = "$(filename)_TailVolumeCoefficient_$(Tail_Volume_Coefficient[i])"
        Tail_Volume_Coefficient[i], DCmg_Dalpha[:, i], Cmg[:, i] = Tail_Analysis(chord_length, wingspan, HS_chord1[i], HS_span1[i], VS_chord, VS_span, tail_distance, wing_distance, n, alpha_range, "$(filename1[i]).csv")
        plot1 = plot!(alpha_range, DCmg_Dalpha[:, i], label= "Tail Volume Coefficient: $(Tail_Volume_Coefficient[i])", xlabel="Alpha(degrees)", ylabel="Stability Derivative")
        savefig(plot1, "$(filename).png")
    end
    plot()
    for i = 1:length(Tail_Volume_Scalars) #plot trim stability for each tail volume scalar
        plot2 = plot!(alpha_range, Cmg[:, i], label= "Tail Volume Coefficient: $(Tail_Volume_Coefficient[i])", xlabel = "Alpha (degrees)", ylabel = "Cmg")
        savefig(plot2, "$(filename)_TrimStability.png")
        CSV_Header = ["Alpha_Range" "Cmg"]
        CSV_Array = Array{Float64, 2}(undef, length(alpha_range), 2)
        CSV_Array[:, 1] = alpha_range
        CSV_Array[:, 2] = Cmg[:, i]
        WriteFile(CSV_Array, "$(filename)_TVC_$(Tail_Volume_Coefficient[i])TrimStability.csv", CSV_Header)
    end


end

function WriteFile(CSVArray, filename, CSVHeader) #write files to CSV
    open(filename, "w") do io
        writedlm(io, CSVHeader, ',')
        writedlm(io, CSVArray, ',')
    end
end

function Derivative_calc(x, y) #inpute the x values and corresponding y values of the function you want to take the first derivative of
    y_prime = Array{Float64, 1}(undef, length(x))
    for i = 1:length(x)
        if i < 3
            y_prime[i] = (-y[i + 2] + 4*y[i+1] - 3*y[i])/(2*(x[i+1]-x[i])) #calculates the forward derivative
        elseif i > 2 && i < length(x) - 1
            y_prime[i] = (-y[i + 2] + 8*y[i+1]-8*y[i-1]+y[i-2])/(12*(x[i]-x[i-1])) #calculates central derivative

        else
            y_prime[i] = (3*y[i]-4*y[i-1]+y[i-2])/(2*(x[i]-x[i-1])) #calculates backwards derivative
        end
    end
    return y_prime
end

function VLM_AlphaRange(Wingspan, chord_length, alpha_range) #This function performs a VLM analysis on a range of alpha values assuming rectangular wingshape
    CFx = Array{Float64, 1}(undef, length(alpha_range))
    CFy = Array{Float64, 1}(undef, length(alpha_range))
    CFz = Array{Float64, 1}(undef, length(alpha_range))
    Cmx = Array{Float64, 1}(undef, length(alpha_range))
    Cmy = Array{Float64, 1}(undef, length(alpha_range))
    Cmz = Array{Float64, 1}(undef, length(alpha_range))
    CDiff = Array{Float64, 1}(undef, length(alpha_range))
    for i = 1:length(alpha_range)
        CFx[i], CFy[i], CFz[i], Cmx[i], Cmy[i], Cmz[i], CDiff[i] = VLM(Wingspan*chord_length, chord_length, Wingspan,alpha_range[i]*pi/180 ,0.0, [0.0; 0.0; 0.0], 1) #Assuming beta and omega are 0 and it also converts alpha_range into radians
    end
    return CFx, CFy, CFz, Cmx, Cmy, Cmz, CDiff #returns the coefficients for force and moment in cartesian coordinates
end

function Pitch_Trim_Stability(Cmw, Clw, Cmt, Clt, wing_distance, tail_distance, alpha_range, Tail_Volume_Coefficient, n, wing_chord) #This function checks trim stability at the given angles of attack. It does this
    #by computing the Moment Coefficient about the center of gravity.
    #It outputs a graph/CSV containing the moment coefficient about the center of gravity.
    #This utilizes equation 5.23 in the flight vehicle design textbook
    #This assumes there is a tail
    #inputs are Moment and lift coefficients for the wing and tail. Also the distance from the aerodynamic center of the tail and wing to the center of gravity.
    #Other inputs are the tail volume coefficient, wing efficiency, and wing chord length
    Cmg = Array{Float64, 1}(undef, length(alpha_range))
    for i = 1:length(Cmg)
        Cmg[i] = Cmw[i] + Cmt[i] - (wing_distance/wing_chord)*Clw[i] - n*Tail_Volume_Coefficient*Clt[i]*tail_distance #This is just equation 5.23 in the flight vehicle design textbook
    end
    return Cmg
end

function VLM(Sref, cref, bref, alpha, beta, omega, Vinf_in) #Performs a Vortex lattice analysis
    xle = [0.0, 0.0] #first number is the position of the leading edge closest to the fuselage in the chordwise direction ,2nd is the same thing but with the leading edge at the wingtip
    yle = [0.0, bref/2] #spanwise direction
    zle = [0.0, 0.0] #vertical direction, use this for dihedral
    chord = [bref/2, bref/2] #first number is chord at the fueslage, the next is the chord at the wingtip.
    theta = [0.0, 0.0] #this is twist (rotation about y-axis) at the fueslage and wingtip respectively.
    phi = [0.0, 0.0] #This is rotation about the x-axis.
    #fc = fill((xc) -> 0, 2) # camberline function for each section, I don't think I need this

    Panels_span = 30
    Panels_chord = 15
    Spacing_type_span = Cosine()
    Spacing_type_chord = Uniform()
    Rref = [0.0,0.0,0.0]
    Vinf = 1
    ref = Reference(Sref, cref, bref, Rref, Vinf)
    fs = Freestream(Vinf, alpha, beta, omega) #Define freestream Parameters

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

function RangeMaker(range, stepsize) #creates an array of values given the range you want and stepsize
    n = round(Int, abs((range[2] - range[1])/stepsize))
    x = Array{Float64, 1}(undef, n)
    for i = 1:n
        if i < 2
            x[1] = range[1]
        elseif i > 1 && i < n
            x[i] = x[i - 1] + stepsize

        else
            x[i] = range[2]
        end
    end
    return x
end

#Tail_Analysis(4, 10, 2, 4, 3, 4, 7, 2, 0.85, RangeMaker([-1 1], 0.01), "Assignment2\\TailTest.csv") #This was for testing the Tail_Analysis funtion
#This is my call that I've been using for research

chord_length = 4
wingspan = 10
HS_chord = 0.5
HS_span = 2
VS_chord = 3 
VS_span = 4
tail_distance = 10
wing_distance = 0
n = 1
alpha_range  = RangeMaker([-1 1], 0.05)
filename = "Assignment2\\Tail_Volume_Research_Docs\\Test10"
Tail_Volume_Scalars = [1 2 3 4 5]
Tail_Volume_Coefficient_Repeater(chord_length, wingspan, HS_chord, HS_span, VS_chord, VS_span, tail_distance, wing_distance, n ,alpha_range, filename, Tail_Volume_Scalars)

println("") #I add this here so it won't print anything unless I need it to