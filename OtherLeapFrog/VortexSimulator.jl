#initialize variables
#Determine starting positions
#build a function that computes the cross product to find induced velocity
using Plots
using LinearAlgebra

@profview begin

function crossproduct(x,y) #cross product function
    z1 = (x[2]*y[3] - y[2]*0)
    z2 = -(x[1]*y[3] - y[1]*0)
    z3 = (x[1]*0 - x[2]*0)
    return z1, z2, z3
end

function RadiusFinder(x,y)
    r = Array{Float16,2}(undef,3,1)
    r[1] = (y[1] - x[1])/2
    r[2] = (y[2] - x[2])/2
    r[3] = 0
    return r
end

global k1 = Array{Float16, 2}(undef, 3, 1)
global k2 = Array{Float16, 2}(undef, 3, 1)
global deltat = 0.01
totalt = 10
n = 4000
global P1 = Array{Float16 , 2}(undef, 3, n) #initialize Point 1
global P2 = Array{Float16 , 2}(undef, 3, n) #initialize Point 2
global P3 = Array{Float16 , 2}(undef, 3, n) #initialize Point 3
global P4 = Array{Float16 , 2}(undef, 3, n) #initialize Point 4

global V1 = Array{Float16 , 2}(undef, 3, n) #initialize Point 1
global V2 = Array{Float16 , 2}(undef, 3, n) #initialize Point 2
global V3 = Array{Float16 , 2}(undef, 3, n) #initialize Point 3
global V4 = Array{Float16 , 2}(undef, 3, n) #initialize Point 4
#define initial conditions

P1[1,1] = 0
P1[2,1] = -0.5
P1[3,1] = 0
P2[1,1] = 0
P2[2,1] = 0.5
P2[3,1] = 0
P3[1,1] = 1
P3[2,1] = 0.5
P3[3,1] = 0
P4[1,1] = 1
P4[2,1] = -0.5
P4[3,1] = 0

V1[1,1] = 0
V1[2,1] = 0
V1[3,1] = 0
V2[1,1] = 0
V2[2,1] = 0
V2[3,1] = 0
V3[1,1] = 0
V3[2,1] = 0
V3[3,1] = 0
V4[1,1] = 0
V4[2,1] = 0
V4[3,1] = 0

global radius1 = Array{Float16,1}(undef,3) #radius from point 1 to point 4
global radius2 = Array{Float16,1}(undef,3) #radius from point 2 to point 3
global radius3 = Array{Float16,1}(undef,3) #radius from point 3 to point 2
global radius4 = Array{Float16,1}(undef,3) #radius from point 4 to point 1

#find the values of the radius, then compute cross product with speed vector
for i = 1:n-1
    global radius1 = RadiusFinder(P1[:,i],P4[:,i])
    global radius2 = RadiusFinder(P2[:,i],P3[:,i])
    global radius3 = RadiusFinder(P3[:,i],P2[:,i])
    global radius4 = RadiusFinder(P4[:,i],P1[:,i])
    R1 = sqrt((2*radius1[1,1])^2 + (2*radius1[2,1])^2)
    R2 = sqrt((2*radius2[1,1])^2 + (2*radius2[2,1])^2)
    global k1 = [0,0, 1/(2*3.141592*R1)]
    global k2 = [0,0, 1/(2*3.141592*R2)]
    V1[1,i + 1],V1[2,i + 1],V1[3,i + 1] = crossproduct(radius4, k1)
    V2[1,i + 1],V2[2,i + 1],V2[3,i + 1] = crossproduct(radius3, k2)
    V3[1,i + 1],V3[2,i + 1],V3[3,i + 1] = crossproduct(radius2, k2)
    V4[1,i + 1],V4[2,i + 1],V4[3,i + 1] = crossproduct(radius1, k1)
    P1[1,i + 1] = P1[1,i] + V1[1,i + 1]*deltat
    P1[2,i + 1] = P1[2,i] + V1[2,i + 1]*deltat
    P1[3,i + 1] = P1[3,i] + V1[3,i + 1]*deltat
    P2[1,i + 1] = P2[1,i] + V2[1,i + 1]*deltat
    P2[2,i + 1] = P2[2,i] + V2[2,i + 1]*deltat
    P2[3,i + 1] = P2[3,i] + V2[3,i + 1]*deltat
    P3[1,i + 1] = P3[1,i] + V3[1,i + 1]*deltat
    P3[2,i + 1] = P3[2,i] + V3[2,i + 1]*deltat
    P3[3,i + 1] = P3[3,i] + V3[3,i + 1]*deltat
    P4[1,i + 1] = P4[1,i] + V4[1,i + 1]*deltat
    P4[2,i + 1] = P4[2,i] + V4[2,i + 1]*deltat
    P4[3,i + 1] = P4[3,i] + V4[3,i + 1]*deltat
end
#=
println("Vecotrs are...")
println(V1[:,2])
println(V2[:,2])
println(V3[:,2])
println(V4[:,2])
=#
staticplot = plot(P1[1,:], P1[2,:])
staticplot = plot!(P2[1,:], P2[2,:])
staticplot = plot!(P3[1,:], P3[2,:])
staticplot = plot!(P4[1,:], P4[2,:])
savefig(staticplot, "Figure1")

end