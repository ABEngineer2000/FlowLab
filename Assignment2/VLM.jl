using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice

#make sure all these numbers aren't integers.
#define varibles
xle = [0.0, 0.0] #first number is the position of the leading edge closest to the fuselage in the chordwise direction ,2nd is the same thing but with the leading edge at the wingtip
yle = [0.0, 10] #spanwise direction
zle = [0.0, 0.0] #vertical direction, use this for dihedral
chord = [2.0, 2.0] #first number is chord at the fueslage, the next is the chord at the wingtip.
theta = [0.0, 0.0] #this is twist (rotation about y-axis) at the fueslage and wingtip respectively.
phi = [0.0, 0.0] #This is rotation about the x-axis.
fc = fill((xc) -> 0, 2) # camberline function for each section

#Lattice spacing
Panels_span = 15
Panels_chord = 4
Spacing_type_span = Uniform()
Spacing_type_chord = Uniform()

#specify Referrence Parameters
Sref = 20.0 #reference surface area
cref = 2.0 #reference chord length
bref = 20.0 #reference span
Rref = [0, 0, 0] #referrence location for all rotations/moments
Vinf = 1.0 #free stream velocity
alpha = 1*pi/180 #pitch angle of attack in radians!
beta = 0.0 # yaw angle of attack
omega = [0.0; 0.0; 0.0] #roll angle of attack

ref = Reference(Sref, cref, bref, Rref, Vinf) #compile reference Parameters
fs = Freestream(Vinf, alpha, beta, omega) #Define freestream Parameters

#create the surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, Panels_span, Panels_chord, spacing_s = Spacing_type_span, spacing_c = Spacing_type_chord)
surfaces = [surface]
#perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric = true)
CF, CM = body_forces(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

println(" ")
println(CF)
println(CM)