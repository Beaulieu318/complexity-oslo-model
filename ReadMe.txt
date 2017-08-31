The code was first written in python then converted to c++. 
I use shared object (dll in Windows) for the C++ code.
These will not work on a MAC but work on windows computers.

To run the code:
Change the working to directory: cd C:\\User\\...

Note: N = number of drives, L = system size, p = probability
      Although the report uses 10^7 drives, this can take up to 10 minutes
      Using 10^6 drives usually takes 1 minutes to get all the data and load all the graphs.

OsloModelPython(N, L, p) #runs the animiation

To get the data, enter this into the terminal:


Plots, heights, avalanches, drops = MutliOsloModelCPP(N=1000000, Lmin=3, Lmax=9, p=0.5)


The follow code can be copied and pasted into the terminal to produce all the relavent graphs for the heights, avalanche-sizes, and drop-sizes.
Please run the code in order.
The specific figures are labelled but should be explicit.

For the heights:

s_heightanalysis = sha.SavedHeightsAnalysis(data = heights, Plots=Plots)
s_heightanalysis.Heights(MovingAverage = False, CalculateCriticals = False, t0=100000, T=600000)
s_heightanalysis.Gradients(points_excluded=0) #this needs to be run
s_heightanalysis.HeightsPlot() #figure1
s_heightanalysis.LogHeightGradientPlot() #figure1
s_heightanalysis.CrossOverTimePlot()
s_heightanalysis.HeightOfPilePlot()
s_heightanalysis.GradientOfPilePlot() #figure3
s_heightanalysis.SDHeightOfPilePlot(alpha=0.2323) #figure4
s_heightanalysis.CorrectionToScalingHeightPlot(a0=1.8, optimise=True)
s_heightanalysis.CorrectionToScalingSDHeightPlot(b0=1.8, alpha=3, optimise=True)
s_heightanalysis.CollapsedHeightsPlot(D=None) #figure2
s_heightanalysis.ProbabilitiesPlot()
s_heightanalysis.CollapsedProbabilitiesPlot() #figure4

For the avalanche-sizes:

s_avalancheanalysis = saa.SavedAvalancheAnalysis(data = avalanches, Plots = Plots)
s_avalancheanalysis.Avalanches(MovingAverage = False, T=700000)
s_avalancheanalysis.DandTauEstimate(a=1.3)
s_avalancheanalysis.RawLogBinAndTauS(a=1.4) #figure5
s_avalancheanalysis.ProbabilitiesPlot(a=1.4) #figure6
s_avalancheanalysis.MomentAnalysis()
s_avalancheanalysis.MomentAnalysisPlot() #figure8
s_avalancheanalysis.DMomentAnalysisPlot() #figure8
s_avalancheanalysis.CollapsedProbabilitiesPlot(a=1.5, D=None, tau_s=None) #figure7

For the drop-sizes:

s_dropanalysis = sda.SavedDropAnalysis(data = drops, Plots = Plots)
s_dropanalysis.Drops(MovingAverage = False)
s_dropanalysis.RawLogBinAndTauS(a=1.6)
s_dropanalysis.ProbabilitiesPlot(a=1.6)
s_dropanalysis.CollapsedProbabilitiesLPlot(a=1.6, alpha=1.2, D=1.25) #figure9
s_dropanalysis.CollapsedProbabilitiesDropPlot(a=1.6, tau_d=1.0, D=1.25) #figure9