using Xfoil, Plots, Printf, LinearAlgebra, DelimitedFiles

function CombinePlots(CSVArrayList, CorrespondingLabelTitles, ExperimentType) #input the list of CSV files you want to plot and the corresponding label titles for each CSV File
    #as well as the Experiment type you are doing aka camber or thickness
    c_l = plot()
    for i = 1:length(CSVArrayList) #for leap for each seperate plot
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_l = Array_plot[2:n, 2]
        lift_to_drag_ratio = Array_plot[2:n, 5]

        plot_lift = plot!(Alpha, c_l, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", grid = false)
        savefig(plot_lift, "Assignment1\\Latex\\$(ExperimentType)Lift_plot.png")
    end
    c_d = plot()
    for i = 1:length(CSVArrayList)
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_d = Array_plot[2:n, 3]

        plot_drag = plot!(Alpha, c_d, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient", grid = false)
        savefig(plot_drag, "Assignment1\\Latex\\$(ExperimentType)Drag_plot.png")
    end
    c_m = plot()
    for i = 1:length(CSVArrayList)
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        c_m = Array_plot[2:n, 4]

        plot_moment = plot!(Alpha, c_m, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient", grid = false)
        savefig(plot_moment, "Assignment1\\Latex\\$(ExperimentType)Moment_plot.png")

    end
    lift_to_drag_ratio = plot()
    for i = 1:length(CSVArrayList)
        Array_plot = readdlm(CSVArrayList[i], ',') #read in csv files
        n = size(Array_plot, 1) - 1 #since the first row will be a string I create an n variable that is the number of rows besides the first
        Alpha = Array_plot[2:n, 1]
        lift_to_drag_ratio = Array_plot[2:n, 5]

        plot_moment = plot!(Alpha, lift_to_drag_ratio, label= CorrespondingLabelTitles[i], xlabel="Angle of Attack (degrees)", ylabel="Lift To Drag Ratio", grid = false)
        savefig(plot_moment, "Assignment1\\Latex\\$(ExperimentType)LiftVsDrag_plot.png")

    end



end

labeltitles = ["2% Camber", "4% Camber", "5% Camber", "9.5% Camber"]
ExperimentType = "Camber"
CSVList = ["Assignment1\\CamberedAirfoils\\2PercentCamber_NACA2412Camber_Percent_2.csv", "Assignment1\\CamberedAirfoils\\4PercentCamber_NACA2412Camber_Percent_4.csv", "Assignment1\\CamberedAirfoils\\5PercentCamber_NACA2412Camber_Percent_5.csv", "Assignment1\\CamberedAirfoils\\9.5PercentCamber_NACA2412Camber_Percent_9.5.csv"]
CombinePlots(CSVList, labeltitles, ExperimentType)