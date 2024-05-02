using Plots #using the plots package, you have to pre install any packages beforehand

function Practice(x,y) #function practice, basically inside the parenthesis is the input of the function, return will return the output
    z = x + y*3
    o = x - 3*y
    return z, o
end

z, o = Practice(3,4)
println(z)
println(Char(' '))
Array1 = readdlm("C:\\Users\\josep\\OneDrive\\Desktop\\FlowLab\\Airfoils\\NACA 16-006.txt", Float16)
x = Array{Float16, 2}(undef, size(Array1, 1), 1)
y = Array{Float16, 2}(undef, size(Array1, 1), 1)
println(size(Array1, 1))
println(length(x))

for i = 1:size(Array1, 1)
    x[i] = Array1[i, 1]
    y[i] = Array1[i, 2]
end
println(x)
println(y)

function Camber(x,c,p)
    if x <= p
        z = c*(2*p*x - x^2)/(p^2)
    else
        z = c*(1-2*p + 2*p*x-x^2)/((1-p)^2)
    end
    return z
end
n = 100 #number of points
x = Array{Float16, 1}(undef, n) #initialize x vector in this case were going from 0 to 1
i = 2
x[1] = 0
for i = 2:n
    x[i] = x[i-1] + 0.01
end
yUpper = Array{Float16, 1}(undef, n)
yLower = Array{Float16, 1}(undef, n)
i = 1
p = 0.4
c = 0.02
#println("x = $x")
z = Array{Float16, 1}(undef, n)
for i = 1:n
    z[i] = Camber(x[i],c,p)
end
#println("z = $z")
i = 1
for i = 1:n
    yUpper[i] = (p/0.2)*(0.2696*sqrt(x[i]) - 0.126*x[i] - 0.3516*(x[i])^2 + 0.2843*(x[i])^3 - 0.1015*(x[i])^4)
    yLower[i] = -(p/0.2)*(0.2696*sqrt(x[i]) - 0.126*x[i] - 0.3516*(x[i])^2 + 0.2843*(x[i])^3 - 0.1015*(x[i])^4)
end
legendhelp = -0.5
plot1 = plot(x,yUpper, xlim = (0,0.9), ylim = (-1,1), label = "NACA2412", linecolor = :steelblue) #I had to divide this into two subplots where one was labeled and the other was unlabeled.
plot1 = plot!(x,yLower, xlim = (0,0.9), ylim = (-1,1), label = false, linecolor = :steelblue)

savefig(plot1, "Plot1.png")

array1 = Array{Int,3}(undef,3,10,2)
for k = 1:2
for i = 1:3
    for j = 1:10
        array1[i,j,k] = i + j*5 + k
    end
end
end
#array1

#a, b = Practice(1,2)
#println("$a")
#println("b = $b") #dollar sign will print the variable
