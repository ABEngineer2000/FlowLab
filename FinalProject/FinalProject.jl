using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, DelimitedFiles, SNOW

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

