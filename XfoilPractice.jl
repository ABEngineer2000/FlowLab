using Xfoil, Plots, Printf


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

