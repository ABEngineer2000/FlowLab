using Xfoil, Plots, Printf, LinearAlgebra, DelimitedFiles

#set range of reynolds number and incremental number
#set alpha range
#read in airfoil coordinates
#complete an xfoil analysis for each reynolds number
#save the results to a txt file or csv file

function ReyonldsExperiment(ReynoldsNumber, AlphaRange, DeltaAlpha, n_pan, iterations)
    global re = ReynoldsNumber
    #read in file"
    Array1 = readdlm("Airfoils\\NACA 16-006.txt", Float16)
    x = Array{Float16, 2}(undef, size(Array1, 1), 1)
    y = Array{Float16, 2}(undef, size(Array1, 1), 1)
    #= #troubleshooting
    println(size(Array1, 1))
    println(length(x))
    =#

    for i = 1:size(Array1, 1)
        x[i] = Array1[i, 1]
        y[i] = Array1[i, 2]
    end
    #=
    x, y = open("C:\\Users\\josep\\OneDrive\\Desktop\\FlowLab\\Airfoils\\NACA 16-006.txt", "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end
    =#

    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)

    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane(npan = n_pan)

    #set range of angles of attack
    n = (AlphaRange[1,2] - AlphaRange[1,1])/DeltaAlpha + 1
    n_int = Int(round(n))
    alpha = Array{Float16, 1}(undef, n_int)
    alpha[1] = AlphaRange[1,1]
    for i = 2:n_int
        alpha[i] = alpha[i - 1] + DeltaAlpha
    end

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)

    #Solve Numerically for Lift, Drag, and Moment
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter= iterations, reinit=true)
    end

    # print results - this is for testing
    #=
    println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
    i = 1
    for i = 1:length(alpha)
        @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
    end
    =#

    #All coefficients are put into an array for exporting
    CSVArray = Array{Float16, 2}(undef, length(alpha), 5)
    CSVArray[:, 1] = alpha
    CSVArray[:, 2] = c_l
    CSVArray[:, 3] = c_d
    CSVArray[:, 4] = c_m
    CSVArray[:, 5] = converged
    CSVHeader = Array{String, 2}(undef, 1, 5)
    CSVHeader[1,:] = ["alpha";"c_l"; "c_d"; "c_m"; "converged"]
     # troubleshooting
     #=
    println(CSVArray[1,:])
    println(CSVArray[2,:])
    println(CSVArray[3,:])
    println(CSVArray[4,:])
    =#

    return alpha, c_l, c_d, c_m, converged, CSVArray, CSVHeader
end

function WriteFile(CSVArray, filename, CSVHeader)
    open(filename, "w") do io
        writedlm(io, CSVHeader, ',')
        writedlm(io, CSVArray, ',')
    end
end

function ReynoldsRepeater(ReynoldsRange, AlphaRange, DeltaAlpha, n_pan, iterations, filename)
    for i = 1:length(ReynoldsRange)
        alpha, c_l, c_d, c_m, converged, CSVArray, CSVHeader = ReyonldsExperiment(ReynoldsRange[i], AlphaRange, DeltaAlpha, n_pan, iterations)
        filename_i = "$(filename)Reynolds_Number$(ReynoldsRange[i]).csv"
        WriteFile(CSVArray, filename_i, CSVHeader)   
    end
end

function plotter(CSVFile, filename)
    Array_plot = readdlm(CSVFile, ',')
    n = size(Array_plot, 1) - 1
    Header = Array_plot[1, :]
    Alpha = Array_plot[2:n, 1]
    c_l = Array_plot[2:n, 2]
    c_d = Array_plot[2:n, 3]
    c_m = Array_plot[2:n, 4]
    c_lplot = plot(Alpha, c_l, label=false, xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient")
    c_dplot = plot(Alpha, c_d, label=false, xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient")
    c_mplot = plot(Alpha, c_m, label=false, xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient")
    savefig(c_lplot, "$(filename)_lift_plot.png")
    savefig(c_dplot,"$(filename)_Drag_plot.png")
    savefig(c_mplot,"$(filename)_Moment_plot.png")
end

plotter("Assignment1\\NACA16-006Reynolds_Number10000.0.csv", "Assignment1\\NACA16-006Reynolds_Number10000.0")

#= 
A = [1 2 3 4; 5 6 7 8; 9 10 11 12]
println(A[:, 2:4 ])
=#

#=
ReynoldsRepeater([1000 10000 1e5 1e6 1e7], [0 3], 0.5, 100, 100, "Assignment1\\NACA16-006")
=#
#call function here
#=
alpha, c_l, c_d, c_m, converged, CSVArray, CSVHeader = ReyonldsExperiment(1e5, [0 3], 0.5, 100, 100)
WriteFile(CSVArray, "Assignment1\\Test1.csv", CSVHeader)
println("")
=#
println("")