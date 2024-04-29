using Xfoil, Plots, Printf, LinearAlgebra


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
rangealpha = [0 14]
deltaalpha = 1
n = (rangealpha[1,2] - rangealpha[1,1])/deltaalpha + 1
n_int = Int(round(n))
alpha = Array{Float16, 1}(undef, n_int)
alpha[1] = rangealpha[1,1]
for i = 2:n_int
    alpha[i] = alpha[i - 1] + deltaalpha
end

re = 1e5 # Reynolds number

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
for i = 1:n_a
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
end

plot1 = plot(alpha, c_l, label = false, xlabel = "Alpha", ylabel = "Coefficient of Lift")
savefig(plot1, "Coefficient of Lift Graph.png")