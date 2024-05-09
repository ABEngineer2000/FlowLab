using Xfoil, Plots, Printf, LinearAlgebra, DelimitedFiles

function CombinePlots(CSVArrayList, CorrespondingLabelTitles, ExperimentType) #input the list of CSV files you want to plot and the corresponding label titles for each CSV File
    #as well as the Experiment type you are doing aka camber or thickness
    plot() #this line is really important, it resets the plot.
    for i = 1:length(CSVArrayList) #for leap for each seperate plot
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_l = Array_plot[2:n, 2]
        lift_to_drag_ratio = Array_plot[2:n, 5]

        plot_lift = plot!(Alpha, c_l, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", grid = false)
        savefig(plot_lift, "Assignment1\\Latex\\$(ExperimentType)Lift_plot.png") #I kept this line in the for loop even though it is inefficient. It is kept here because if I put it outside the code won't be able to access the variable. I'd rather leave it in than deal with the hassle of making it a global varible.
    end
    plot()
    for i = 1:length(CSVArrayList)
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_d = Array_plot[2:n, 3]

        plot_drag = plot!(Alpha, c_d, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient", grid = false)
        savefig(plot_drag, "Assignment1\\Latex\\$(ExperimentType)Drag_plot.png")
    end
    plot()
    for i = 1:length(CSVArrayList)
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_m = Array_plot[2:n, 4]

        plot_moment = plot!(Alpha, c_m, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient", grid = false)
        savefig(plot_moment, "Assignment1\\Latex\\$(ExperimentType)Moment_plot.png")

    end
    plot()
    for i = 1:length(CSVArrayList)
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        lift_to_drag_ratio = Array_plot[2:n, 5]

        plot_moment = plot!(Alpha, lift_to_drag_ratio, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Lift To Drag Ratio", grid = false)
        savefig(plot_moment, "Assignment1\\Latex\\$(ExperimentType)LiftVsDrag_plot.png")

    end



end

labeltitles = ["Generated Data,", "Airfoil Tools Data"]
ExperimentType = "Angle_Of_Attack_Experiment"
CSVList = ["Assignment1\\NACA16-006\\NACA16-006Reynolds_Number10000.0.csv",""]
#CombinePlots(CSVList, labeltitles, ExperimentType)

#= #this is for plotting something very specific, probably won't need to use again.
Array_plot = readdlm("Assignment1\\NACA16-006\\NACA16-006Reynolds_Number10000.0.csv", ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_l = Array_plot[2:n, 2]
        AOE_Lift = plot(Alpha, c_l, label= "Generated Data", xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", grid = false)
Array_plot2 = readdlm("Assignment1\\ComparisonData\\NACA16-006_LIft.csv",',')
    Alpha = Array_plot2[:, 1]
    c_l = Array_plot2[:, 2]
    AOE_Lift = plot!(Alpha, c_l, label= "Airfoil Tools Data", xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", grid = false)
    savefig(AOE_Lift, "Assignment1\\Latex\\AOE_Lift_plot.png")

    plot()
Array_plot = readdlm("Assignment1\\NACA16-006\\NACA16-006Reynolds_Number10000.0.csv", ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_d = Array_plot[2:n, 3]
        AOE_drag = plot(Alpha, c_d, label= "Generated Data", xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient", grid = false)
Array_plot = readdlm("Assignment1\\ComparisonData\\NACA16-006_Drag.csv",',')
    Alpha = Array_plot[:, 1]
    c_d = Array_plot[:, 2]
    AOE_drag = plot!(Alpha, c_d, label= "Airfoil Tools Data", xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient", grid = false)
    savefig(AOE_drag, "Assignment1\\Latex\\AOE_drag_plot.png")

    
    plot()
Array_plot = readdlm("Assignment1\\NACA16-006\\NACA16-006Reynolds_Number10000.0.csv", ',') #read in csv files
    n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
    Alpha = Array_plot[2:n, 1]
    c_m = Array_plot[2:n, 4]
    AOE_moment = plot(Alpha, c_m, label= "Generated Data", xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient", grid = false)
Array_plot2 = readdlm("Assignment1\\ComparisonData\\NACA16-006_Moment.csv",',')
    Alpha = Array_plot2[:, 1]
    c_m = Array_plot2[:, 2]
    AOE_moment = plot!(Alpha, c_m, label= "Airfoil Tools Data", xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient", grid = false)
    savefig(AOE_moment, "Assignment1\\Latex\\AOE_Moment_plot.png")
    =#