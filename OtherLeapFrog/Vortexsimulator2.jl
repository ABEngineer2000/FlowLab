using Plots
using LinearAlgebra

#initialize variables
#define initial conditions
#compute distance between each point for each point
#create velocity vector from those distances
#compute change in position for each point

#=
Cross product a = c x d, a = cross_op(c)*d
https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
=#

#= Practice with matrices
A = [1 2; 3 4]
B = [2 0; 0 2]
C = transpose(A)*B
=#

function cross_product(x,y) #cross operator matrix, x is the left vector in the cross product operation
    x_op = [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]
    z = x_op*y
    return z
end

function DistanceFinder(x,y)
    r = [(y[1] - x[1]); (y[2] - x[2]); 0]
    return r
end

function Mag(x)
    mag = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    return mag
end

#=
#testing the cross product function
c = [1;2;3]
d = [1;0;4]

a = cross_product(c,d)
=#
#define variables and initial conditions\
deltat = 0.01
n = 4000
global P1 = Array{Float16 , 2}(undef, 3, n) #initialize Point 1
global P2 = Array{Float16 , 2}(undef, 3, n) #initialize Point 2
global P3 = Array{Float16 , 2}(undef, 3, n) #initialize Point 3
global P4 = Array{Float16 , 2}(undef, 3, n) #initialize Point 4

P1[:,1] = [0;-0.5;0]
P2[:,1] = [0;0.5;0]
P3[:,1] = [1;0.5;0]
P4[:,1] = [1;-0.5;0]

global k = [0;0;1]
#main for loop to calculate velocity and position
for i = 1:n-1
    local v1 = [0;0;0]
    local v2 = [0;0;0]
    local v3 = [0;0;0]
    local v4 = [0;0;0]
    v1 = cross_product(-k,DistanceFinder(P1[:,i],P2[:,i]))/(2*3.1415926*Mag(DistanceFinder(P1[:,i],P2[:,i]))^2) + cross_product(-k,DistanceFinder(P1[:,i],P3[:,i]))/(2*3.1415926*Mag(DistanceFinder(P1[:,i],P3[:,i]))^2) + cross_product(k,DistanceFinder(P1[:,i],P4[:,i]))/(2*3.1415926*Mag(DistanceFinder(P1[:,i],P4[:,i]))^2)
    v2 = cross_product(k,DistanceFinder(P2[:,i],P1[:,i]))/(2*3.1415926*Mag(DistanceFinder(P2[:,i],P1[:,i]))^2) + cross_product(-k,DistanceFinder(P2[:,i],P3[:,i]))/(2*3.1415926*Mag(DistanceFinder(P2[:,i],P3[:,i]))^2) + cross_product(k,DistanceFinder(P2[:,i],P4[:,i]))/(2*3.1415926*Mag(DistanceFinder(P2[:,i],P4[:,i]))^2)
    v3 = cross_product(-k,DistanceFinder(P3[:,i],P2[:,i]))/(2*3.1415926*Mag(DistanceFinder(P3[:,i],P2[:,i]))^2) + cross_product(k,DistanceFinder(P3[:,i],P1[:,i]))/(2*3.1415926*Mag(DistanceFinder(P3[:,i],P1[:,i]))^2) + cross_product(k,DistanceFinder(P3[:,i],P4[:,i]))/(2*3.1415926*Mag(DistanceFinder(P3[:,i],P4[:,i]))^2)
    v4 = cross_product(-k,DistanceFinder(P4[:,i],P2[:,i]))/(2*3.1415926*Mag(DistanceFinder(P4[:,i],P2[:,i]))^2) + cross_product(-k,DistanceFinder(P4[:,i],P3[:,i]))/(2*3.1415926*Mag(DistanceFinder(P4[:,i],P3[:,i]))^2) + cross_product(k,DistanceFinder(P4[:,i],P1[:,i]))/(2*3.1415926*Mag(DistanceFinder(P4[:,i],P1[:,i]))^2)
    #=
    debugging
    println(v1)
    println(v2)
    println(v3)
    println(v4)
    =#
    
    P1[:, (i + 1)] = P1[:,i] + v1*deltat
    P2[:, (i + 1)] = P2[:,i] + v2*deltat
    P3[:, (i + 1)] = P3[:,i] + v3*deltat
    P4[:, (i + 1)] = P4[:,i] + v4*deltat
end
x = [1,2,3]
y= [4,5,6]
cross_product(x,y)
@time cross_product(x, y)
#=
anim = Animation()
staticplot = scatter( ;xlim = (0,10), ylim = (-1,1), legend = false)
anim = @animate for j = 1:n
scatter!((P1[1,j],P1[2,j]))
scatter!((P2[1,j],P2[2,j]))
scatter!((P3[1,j],P3[2,j]))
scatter!((P4[1,j],P4[2,j]))
frame(anim, staticplot)
end
=#

#=
savefig(staticplot, "Figure2")
gif(anim, "anim_fps150.gif", fps = 150)
=#