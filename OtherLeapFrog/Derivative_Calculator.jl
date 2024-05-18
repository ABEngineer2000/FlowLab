using Printf, LinearAlgebra, DelimitedFiles

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

function e_x(x_Range, x_step) #inputs e^x for every x in x_Range with a given x_Step size
    n = round(Int, abs((x_Range[2] - x_Range[1])/x_step))
    println(n)
    x = Array{Float64, 1}(undef, n)
    e_x = Array{Float64, 1}(undef, length(x))
    for i = 1:length(x)
        if i < 2
            x[i] = 1
        else
        x[i] = x[i-1] + x_step
        e_x[i] = exp(x[i])
        end
    end
    return x, e_x
end

 x, ex = e_x([1 10], 0.01)
 y_prime = Derivative_calc(x, ex)

 println(x[900])
 println(y_prime[900])

println("")
