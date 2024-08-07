using Xfoil, Plots, Printf, LinearAlgebra

function Xfoil_function(alpharange, dalpha, re)
# read airfoil coordinates from a file
x, y = open("Airfoil1.txt", "r") do f
    x = Float64[]
    y = Float64[]
    for line in eachline(f)
        entries = split(chomp(line))
        push!(x, parse(Float64, entries[1]))
        push!(y, parse(Float64, entries[2]))
    end
    x, y
end



# load airfoil coordinates into XFOIL
Xfoil.set_coordinates(x,y)

# plot the airfoil geometry
scatter(x, y, label="", framestyle=:none, aspect_ratio=1.0, show=true)

# repanel using XFOIL's `PANE` command
xr, yr = Xfoil.pane(npan = 100)


#set range of angles of attack
rangealpha = alpharange
deltaalpha = dalpha
n = (rangealpha[1,2] - rangealpha[1,1])/deltaalpha + 1
n_int = Int(round(n))
alpha = Array{Float16, 1}(undef, n_int)
alpha[1] = rangealpha[1,1]
for i = 2:n_int
    alpha[i] = alpha[i - 1] + deltaalpha
end

# plot the refined airfoil geometry
scatter(xr, yr, label="", framestyle=:none, aspect_ratio=1.0, show=true)

# initialize outputs
n_a = length(alpha)
c_l = zeros(n_a)
c_d = zeros(n_a)
c_dp = zeros(n_a)
c_m = zeros(n_a)
converged = zeros(Bool, n_a)
i = 1
for i = 1:n_a
    c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
end

# print results
println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
i = 1
for i = 1:length(alpha)
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
end

    return  c_l, c_d, c_dp, c_m, converged, alpha
end

function plotter(alpha, c_l, c_d, c_m)
plot1 = plot(alpha, c_l, label="", xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient")
plot2 = plot(alpha, c_d, label="", xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient")
plot3 = plot(alpha, c_m, label="", xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient")
    return plot1, plot2, plot3
end

alpharange = [-9 14]
dalpha = 1
re = 100000
file = "C:\\Users\\josep\\OneDrive\\Desktop\\FlowLab\\Airfoils\\NACA 16-006.txt"
c_l, c_d, c_dp, c_m, converged, alpha = Xfoil_function(alpharange, dalpha, re)
println("")
