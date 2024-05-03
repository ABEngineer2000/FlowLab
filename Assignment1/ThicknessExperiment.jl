using Xfoil, Plots, Printf, LinearAlgebra, DelimitedFiles

function AirfoilAnalysis(AirfoilDat, AlphaRange, DeltaAlpha, n_pan, iterations)
    
    global re = 1e5
    
    #read in Airfoil Data and input into x and y vectors
    Array1 = readdlm(AirfoilDat, Float16)
    x = Array{Float16, 2}(undef, size(Array1, 1), 1)
    y = Array{Float16, 2}(undef, size(Array1, 1), 1)

    for i = 1:size(Array1, 1)
        x[i] = Array1[i, 1]
        y[i] = Array1[i, 2]
    end
    #load in airfoil data in xfoil and re pan the airfoil
    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)

    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane(npan = n_pan)

    #create an Alpha range
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
    ltd = zeros(n_a)
    converged = zeros(Bool, n_a)

    #Solve Numerically for Lift, Drag, and Moment
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter= iterations, reinit=true)
    end

    #solve for lift_to_drag_ratio
    for i = 1:n_a
        ltd[i] = c_l[i]/c_d[i]
    end

    #All coefficients are put into an array for exporting
    CSVArray = Array{Float16, 2}(undef, length(alpha), 6)
    CSVArray[:, 1] = alpha
    CSVArray[:, 2] = c_l
    CSVArray[:, 3] = c_d
    CSVArray[:, 4] = c_m
    CSVArray[:, 5] = ltd
    CSVArray[:, 6] = converged
    CSVHeader = Array{String, 2}(undef, 1, 6)
    CSVHeader[1,:] = ["alpha";"c_l"; "c_d"; "c_m"; "ltd"; "converged"]


    return alpha, c_l, c_d, c_m, ltd, converged, CSVArray, CSVHeader   
end

function WriteFile(CSVArray, filename, CSVHeader)
    open(filename, "w") do io
        writedlm(io, CSVHeader, ',')
        writedlm(io, CSVArray, ',')
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
    lift_to_drag_ratio = Array_plot[2:n, 5]
    lift_to_drag_ratio_plot = plot(Alpha, lift_to_drag_ratio, label=false, xlabel="Angle of Attack (degrees)", ylabel="Lift To Drag Ratio")
    c_lplot = plot(Alpha, c_l, label=false, xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient")
    c_dplot = plot(Alpha, c_d, label=false, xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient")
    c_mplot = plot(Alpha, c_m, label=false, xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient")
    savefig(c_lplot, "$(filename)_lift_plot.png")
    savefig(c_dplot,"$(filename)_Drag_plot.png")
    savefig(c_mplot,"$(filename)_Moment_plot.png")
    savefig(lift_to_drag_ratio_plot, "$(filename)_LiftToDrag_plot.png")
    
end

#initialize variables will end up changes these for each test
camber = 3
thickness = 30
AlphaRange = [-9 16]
DeltaAlpha = 1
n_pan = 100
iterations = 100
filename = "Assignment1\\ChangingThicknessAirfoils\\$(thickness)PercentThickness_NACA2412"
AirfoilDat = "Assignment1\\ChangingThicknessAirfoils\\$(thickness)PercentThickness_NACA2412.txt"


#call functions here
alpha, c_l, c_d, c_m, ltd, converged, CSVArray, CSVHeader = AirfoilAnalysis(AirfoilDat, AlphaRange, DeltaAlpha, n_pan, iterations)
println("")
filename_i = "$(filename).csv"
WriteFile(CSVArray, filename_i, CSVHeader)
plotter("$(filename).csv", "$(filename)_plot.png")