using Xfoil, Plots, Printf, LinearAlgebra, DelimitedFiles

#set range of reynolds number and incremental number
#set alpha range
#read in airfoil coordinates
#complete an xfoil analysis for each reynolds number
#save the results to a txt file or csv file

function ReyonldsExperiment(ReynoldsRange, AlphaRange, DeltaAlpha, n_pan, iterations)
    global re = 1e5
    #read in file"
    Array1 = readdlm("C:\\Users\\josep\\OneDrive\\Desktop\\FlowLab\\Airfoils\\NACA 16-006.txt", Float16)
    x = Array{Float16, 2}(undef, size(Array1, 1), 1)
    y = Array{Float16, 2}(undef, size(Array1, 1), 1)
    println(size(Array1, 1))
    println(length(x))

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
    CSVArray = Array{Float16, 2}(undef, 4, length(alpha))
    CSVArray[1, :] = c_l
    CSVArray[2, :] = c_d
    CSVArray[3, :] = c_m
    CSVArray[4, :] = converged
     # troubleshooting
     #=
    println(CSVArray[1,:])
    println(CSVArray[2,:])
    println(CSVArray[3,:])
    println(CSVArray[4,:])
    =#

    return alpha, c_l, c_d, c_m, converged, CSVArray
end

function WriteFile(CSVArray)

end

#call function here
ReyonldsExperiment(1, [0 3], 0.5, 100, 100)
println("")