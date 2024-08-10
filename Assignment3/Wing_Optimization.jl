using Plots, Printf, LinearAlgebra, DelimitedFiles, VortexLattice, SNOW

#Plan of Action
#The induced drag from VLM will be the function I am optimizing.
#Define a new function that calls VLM. It will be written in the
#form that SNOW can use
#Write all my constraints using equation 5.1 in the optimization textbook
#Set my f as my induced drag
#set lift force, pitch angle, wingspan, and speed as constraints

#tips to add
#reverse AD for the first line of f then forward AD for g when it comes to taking derivatives.

function Wing_Optimization!(g, x)
    leading_edge_distribution = Array{Float64, 1}(undef, length(x))
    #this for loop keeps all the quarter chords aligned. One of my foced constraints
    for i = 1:length(x)
        if i < 2
            leading_edge_distribution[i] = 0.0
        else
           leading_edge_distribution[i] = x[i - 1] - (x[i-1] - x[i])*0.25 #start at the previous chord position and add whatever the difference between the quarter chord lengths there is.
        end
        #=
           if x[i] > x[i - 1] #I added this in here to force the chords to follow a monotinicity pattern because
            #it was ignoring that in my constraint
            x[i] = x[i] - (x[i] - x[i - 1])*1.05
        end
        =#
    end
    
    #println(x)
    CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area = VLM(leading_edge_distribution, x, 8.0)
    #println(0.5*1.007*wing_area*CFz)
    #g = Array{Float64, 1}(undef, length(x))
    g[1] = -0.5*1.007*wing_area*CFz + 1.7 #minimum lift must be 1.7 newtons. I'm deciding to put bounds on this constraint when I call the function.
    #here's my monotinicity function to decrease the chord as x progresses
    for i = 2:length(x0)
        g[i] = x[i] - x[i - 1]
    end
    #=
    g[2] = x[2] - x[1] + 0.0
    println(x[2] - x[1]) 
    g[3] = x[3] - x[2] + 0.0
    =#
    #=
    for i = 2:length(x)
        g[i] = x[i] - x[i - 1]
    end
    =#
    
    #this is just the lift equation using coefficient of lift. The density of air is 1.007 m3/kg for an alititude of 2000 meters

    #debugging
    #println(0.5*1.007*wing_area*CFz)
    #println(CDiff)
    return CDiff
end

#This function will input an x vector with the corresponding y values and create an nth order polyomial fit
function Poly_Regression(x, y, n)
    A = Array{Float64, 1}(undef, n + 1) #A Matrix with the coefficients of the polynomial
    x_matrix = Array{Float64, 2}(undef, length(x), n + 1)
    for i = 1:n + 1
        for j = 1:length(x)
            x_matrix[j, i] = x[j]^(i - 1)
        end
    end
    A = inv(transpose(x_matrix)*x_matrix)*transpose(x_matrix)*y
    #A = [1 2 3 4]
    return A
end

function Wing_smoother(xopt, deltax) #This function outputs a smoothed wing based on polynomial Poly_Regression
    #plan for this function
    #input chord lengths and deltax
    #creates positional and chord length vecotrs
    #inputs those values as x and y into Poly_Regression which will output the polynomial coefficients
    #creates a new chord length vector which has been smoothed
    #outputs the new chord length vector with the corresponding positional vector

    #create theh position vector x, xopt is the chord length
    x = Array{Float64, 1}(undef, length(xopt))
    for i = 1:length(xopt)
        if i < 2
            x[i] = 0.0
        else
            x[i] = x[i - 1] + deltax
        end
    end

    #get the coefficients for the regression fit
    A = Poly_Regression(x, xopt, 2)

    chord_new = Array{Float64, 1}(undef, length(xopt))
    for i = 1:length(xopt)
        chord_new[i] = A[1] + A[2]*x[i] + A[3]*x[i]^2
    end
    return x, chord_new, A

end

#this function performs a finalized wing analysis given the smoothed wing
function Finalized_Wing_Analysis(A, n)
    #all that it needs for an input is the coefficients of the A Matrix and the number of panels you want
    x = Array{Float64, 1}(undef, n + 1)
    deltax = 4.0 / (n+1)
    for i = 1:n + 1
        if i < 2
            x[i] = A[1] #Evaluating the polynomial at 0
        else
            x[i] = A[1] + (deltax*i)*A[2] + A[3]*(deltax*i)^2 #Evaluting the polynomial at x
        end
    end
    leading_edge_distribution = Array{Float64, 1}(undef, length(x))
    #this for loop keeps all the quarter chords aligned. One of my foced constraints
    for i = 1:length(x)
        if i < 2
            leading_edge_distribution[i] = 0.0
        else
           leading_edge_distribution[i] = x[i - 1] - (x[i-1] - x[i])*0.25 #start at the previous chord position and add whatever the difference between the quarter chord lengths there is.
        end
    end
    CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area, Γ = VLM(leading_edge_distribution, x, 8.0, Γ_output = true)
    lift_force = 0.5*1.007*wing_area*CFz
    return CFz, CDiff, wing_area, lift_force, Γ
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

function VLM(leading_edge_distribution, chord_distribution, span; Γ_output = false) #Performs a Vortex lattice analysis
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

    Panels_span = length(xle)
    Panels_chord = 10
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
    Gamma = (system.Γ)[:] #pull out Γ values for plotting lift distribution.
    CFx, CFy, CFz = CF
    CMx, CMy, CMz = CM
    if Γ_output == false
        return CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area #, grid, surface these I added as outputs for troubleshooting
    else
        return CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area, Gamma
    end
end

function leading_edge_finder(x) #this function finds the leading edge coordinates based on the chord length that was optimized
    leading_edge_distribution = Array{Float64, 1}(undef, length(x))
    for i = 1:length(x)
        if i < 2
            leading_edge_distribution[i] = 0.0
        else
           leading_edge_distribution[i] = leading_edge_distribution[i - 1] - (x[i-1] - x[i])*0.25 #start at the previous chord position and add whatever the difference between the quarter chord lengths there is.
        end
    end
    return leading_edge_distribution
end

function wing_plotter(xopt, leading_edge_distribution, wingspan, filename) #this function plots the wing generated
    y_upper = Array{Float64, 1}(undef, length(xopt))
    y_lower = Array{Float64, 1}(undef, length(xopt))
    x = Array{Float64, 1}(undef, length(xopt))
    deltax = wingspan / ((length(leading_edge_distribution) - 1)*2)
    for i = 1:length(xopt)
        if i < 2
            x[i] = 0.0
            y_upper[i] = leading_edge_distribution[i]
            y_lower[i] = leading_edge_distribution[i] - xopt[i]
        else
            x[i] = x[i - 1] + deltax
            y_upper[i] = leading_edge_distribution[i]
            y_lower[i] = leading_edge_distribution[i] - xopt[i]
        end
    end
    plot() #resets plot
    plot1 = plot(x, y_upper, label = "Leading Edge Coordinate")
    plot1 = plot!(x, y_lower, label = "Trailing Edge Coordinate")
    savefig(plot1, "$(filename)_Wingplot.png")
end

function lift_distribution_plotter(Γ, wingspan, filename)
    plot() #resets plot
    x = Array{Float64, 1}(undef, length(Γ)) #define x values for plotting
    for i = 1:length(Γ)
        if i == 1
            x[i] = -4.0
        else
            x[i] = x[i - 1] + wingspan / (length(Γ))
        end 
    end
    plot1 = plot(x, Γ, xlabel="Wing Position", ylabel="Γ value (scalar multiple of lift)", label="Γ Distribution")
    savefig(plot1, filename)
    plot() #resets plot
end

#defining variables
pitch_angle = 5*pi/180 #set pitch angle to 5 degrees
min_lift = 1.7 #set the minimum lift weight to 1.7 Newtons.
wing_span = 8.0 #set the wingspan to 8 meters
speed = 1.0 #set the speed to 1 m/s
density = 1.007 #This is in kg/m3

#setting bounds to pass into the function
x0 = [6.9, 6.8, 6.75, 6.74, 6.73, 6.72, 6.71, 6.70, 6.69, 6.5]

lx = Array{Float64, 1}(undef, length(x0))
ux = Array{Float64, 1}(undef, length(x0))
#=
lg = Array{Float64, 1}(undef, length(x0) + 1) #the plus 1 is for the extra constraint
ug = Array{Float64, 1}(undef, length(x0) + 1)
=#
# ----- set some options ------
ip_options = Dict(
    "max_iter" => 300,
    "tol" => 1e-6
)
solver = IPOPT(ip_options)
options = Options(;solver)
for i = 1:length(x0)
        lx[i] = 0.1 #minimum chord length is 0.1 meters
        ux[i] = 7.0 #maximum chord length is 7.0 meters
end


ng = length(x0) #one extra constraint for lift, then a constraint for each chord length decreasing from the previous one
lg = -Inf*ones(ng)
ug = 0*ones(ng)

#Here's where I exectue the optimization


xopt, fopt, info = minimize(Wing_Optimization!, x0, ng, lx, ux, lg, ug, options)
println(xopt)
leading_edge_distribution = leading_edge_finder(xopt)
wing_plotter(xopt, leading_edge_distribution, 8.0, "Assignment3\\Test15_BeforeSmoothing_100Iterations")
deltax = (4/(length(xopt) - 1))
x, y, A = Wing_smoother(xopt, deltax)
xnew = leading_edge_finder(y)
wing_plotter(y, xnew, 8.0, "Assignment3\\Test15_100Iterations")
CFz, CDiff, wing_area, lift_force, Γ = Finalized_Wing_Analysis(A, 15)
lift_distribution_plotter(Γ, 8.0, "Assignment3\\Test15_Lift_Distribution.png")
println("Smoothed Wing Numerics")
println("Lift Force = $(lift_force)")
println("CDiff = $(CDiff)")


#end optimization execution

#ng = 2
#xopt, fopt, info = minimize(Wing_Optimization!, x0, ng, lx, ux)
#This is a working snow example using the fx function
#=
x0 = [-0.5; -0.5]
ng = 2
xopt, fopt, info = minimize(fx!, x0, ng)
=#

#=
A = Poly_Regression([1, 2, 3, 4, 5, 6, 7], [1; 4; 8; 15; 16; 17; 19], 2)
println(A)
=#


#CFx, CFy, CFz, CMx, CMy, CMz, CDiff, wing_area = VLM([0.0 0.0], [3.0 3.0], 8)
println("")
#println(CFz)
#println(wing_area)
#println(surface)

#=
Interesting results

No matter how I set my initial conditions, the optimizer always tries to bring the chord length to its max value at the point of 
the wing connected to the fueslage. Test 12 illustrates this

The more rectangular the initial conditions then the shape is more uniform in the center section like in Test 10. You
get more of a parabolic shape


=#