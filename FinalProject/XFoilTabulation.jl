using Plots, Printf, LinearAlgebra, DelimitedFiles, Xfoil

#this function inputs the airfoil txt source file as well as the desire
#file output name and then outputs a tabulate value of lift and drag coefficients for
#the 2d airfoil at specified angles of attack at the specified reynolds number
function Xfoil_Tabulation(Airfoil_source, filename, Reynolds_Number, panels_number, iterations)

    #create the angles of attack. It has a resolution of 0.01 degrees and goes from -60 to 60 degrees
    AOA = Array{Float64, 1}(undef, 12000)
    for i = 1:length(AOA)
        if i < 2
            AOA[i] = -60.0
        else
            AOA[i] = AOA[i - 1] + 0.01
        end
    end

    #gather the x and y airfoil data from the airfoil source
    AirfoilData = readdlm(Airfoil_source, Float64)
    x = Array{Float16, 1}(undef, size(AirfoilData, 1))
    y = Array{Float16, 1}(undef, size(AirfoilData, 1))
    for i = 1:length(x)
        x[i] = AirfoilData[i, 1]
        y[i] = AirfoilData[i, 2]
    end

    #load coordinates into Xfoil
    Xfoil.set_coordinates(x,y)

    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane(npan = panels_number)

    # initialize outputs
    n_a = length(AOA)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)

    #Solve Numerically for Lift, Drag, and Moment
    for i = 1:length(AOA)
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(AOA[i], Reynolds_Number; iter= iterations, reinit=true)
        println(AOA[i])
    end

    #All coefficients are put into an array for exporting
    CSVArray = Array{Float16, 2}(undef, length(AOA), 5)
    CSVArray[:, 1] .= AOA
    CSVArray[:, 2] .= c_l
    CSVArray[:, 3] .= c_d
    CSVArray[:, 4] .= c_m
    CSVArray[:, 5] .= converged
    CSVHeader = Array{String, 2}(undef, 1, 5)
    CSVHeader[1,:] = ["alpha";"c_l"; "c_d"; "c_m"; "converged"]

    #write to file
    WriteFile(CSVArray, filename, CSVHeader)

end

function WriteFile(CSVArray, filename, CSVHeader)
    open(filename, "w") do io
        writedlm(io, CSVHeader, ',')
        writedlm(io, CSVArray, ',')
    end
end

#probably need at least 300 iterations to get convergence
Xfoil_Tabulation("FinalProject\\Airfoils\\NACA_6412.txt", "FinalProject\\Tabulated_Airfoil_Data\\NACA_6412.csv", 1000, 50, 300)
println("Done")

