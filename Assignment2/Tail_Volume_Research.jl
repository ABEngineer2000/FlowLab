using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, FiniteDiff, Ranges

function Tail_Analysis(chord_length, wingspan, HS_chord, HS_span, VS_chord, VS_span, tail_distance, wing_distance, n, alpha_range) #Assuming you have a rectangular wing and tail. HS 
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

        VS_dclz_dalpha = Derivative_calc(alpha_range, vCFz)
        VS_dcmy_dalpha = Derivative_calc(alpha_range, vCMy)

        #compute tail volume coefficient. See page 94 Flight Vehicle Design
        Tail_Volume_Coefficient = (HS_span * HS_chord)/(wingspan * chord_length)
        #compute stability derivatives. Note that this is assuming the cmy moments oppose the moments created by clz
        DCmg_Dalpha = Array{Float64, 1}(undef, length(wing_dclz_dalpha))
        for i = 1:length(wing_dclz_dalpha)
            DCmg_Dalpha[i] = wing_dcmy_dalpha[i] + HS_dcmy_dalpha[i]  - (wing_distance/chord_length)*wing_dclz_dalpha[i] - n*Tail_Volume_Coefficient*HS_dclz_dalpha[i]
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

function plotter(CSV, filename) #This function will plot the results based on a csv file and specified filename


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
        CFx[i], CFy[i], CFz[i], Cmx[i], Cmy[i], Cmz[i], CDiff[i] = VLM(Wingspan*chord_length, chord_length, Wingspan,alpha_range[i],0.0, [0.0; 0.0; 0.0], 1) #Assuming beta and omega are 0
    end
    return CFx, CFy, CFz, Cmx, Cmy, Cmz, CDiff #returns the coefficients for force and moment in cartesian coordinates

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
    Panels_chord = 20
    Spacing_type_span = Cosine()
    Spacing_type_chord = Uniform()
    Rref = [0,0,0]
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

#x = collect(1:0.5:5) #this function tests ranges
#println(x)


println("") #I add this here so it won't print anything unless I need it to