import numpy as np
import ctypes
import OsloModel as om
import SavedHeightsAnalysis as sha
import SavedAvalancheAnalysis as saa
import SavedDropAnalysis as sda
import os

def OsloModelPython(N, L, p):
    """N = drives, L = system size, p = probability"""
    model = om.OsloModel(L)
    model.Iteration(N=N, p=p, Box=True)
    return N, L, model.critical, model.heights, model.avalanchesize, model.dropsize
    
def OsloModelCPP(N, L, p):
    """N = drives, L = system size, p = probability"""
    testlib = ctypes.cdll.LoadLibrary(os.getcwd()+"\OsloModelPy.dll")
    critical = np.zeros(1, np.int)
    system = np.zeros(L+1,np.int)
    rands = np.zeros(L+1, np.int)
    avalancheSize = np.zeros(N, np.int)
    dropSize = np.zeros(N, np.int)
    heights = np.zeros(N, np.int)
    testlib.iterations.restype = None
    testlib.iterations(ctypes.c_int(N), ctypes.c_int(L), ctypes.c_float(p), np.ctypeslib.as_ctypes(system), np.ctypeslib.as_ctypes(rands), np.ctypeslib.as_ctypes(avalancheSize), 
                        np.ctypeslib.as_ctypes(dropSize), np.ctypeslib.as_ctypes(heights), np.ctypeslib.as_ctypes(critical))
    critical = critical[0]
    return N, L, critical, heights, avalancheSize, dropSize
    
def MutliOsloModelCPP(N, Lmin, Lmax, p=0.5):
    """runs the Oslo model and saves it to arrays. These arrays can be used in the analysis of the saved data"""
    L_range = [2**i for i in range(Lmin, Lmax+1)]
    All_Data = []
    Plots = []
    data = []
    for L_num in range(len(L_range)):
        All_Data = OsloModelCPP(N, L_range[L_num], p)
        Plots.append(All_Data[:3])
        data.append(All_Data[3:])
        print "L: ", L_range[L_num]
    Plots = np.array(Plots)
    data = np.rollaxis(np.array(data),1)
    height, avalanche, drop = data
    return Plots, height, avalanche, drop
    
def UncertaintiesSDHOP(R, N, Lmin, Lmax, p=0.5):
    """Can run repeat readings to find the uncertainty on the standard deviation of the heights"""
    SDHOP = np.zeros((R, Lmax-Lmin+1))
    for r in range(R):
        Plots, heights, avalanches, drops = MutliOsloModelCPP(N, Lmin, Lmax, p)
        s_heightanalysis = sha.SavedHeightsAnalysis(data=heights, Plots=Plots)
        s_heightanalysis.Heights(MovingAverage = False, CalculateCriticals = False, t0=N*0.1, T=N-N*0.4)
        SDHOP[r] = np.array(s_heightanalysis.SDHOPs)
    SDHOP = np.rollaxis(SDHOP,1)
    SDHOP_data = np.array([Plots[:,1],np.mean(SDHOP, 1), np.std(SDHOP, 1)]/np.sqrt(R))
    s_heightanalysis.Gradients(points_excluded=4)
    s_heightanalysis.SDHeightOfPilePlot(alpha=0.2323, data=SDHOP_data)
    
def UncertaintiesMoment(R, N, Lmin, Lmax, p=0.5):
    """Can run repeat readings to find the uncertainty on the value of the moments"""
    Moment = np.zeros((R, Lmax-Lmin+1))
    DTau = np.zeros((R, 2))
    for r in range(R):
        Plots, heights, avalanches, drops = MutliOsloModelCPP(N, Lmin, Lmax, p)
        s_avalancheanalysis = saa.SavedAvalancheAnalysis(data=avalanches, Plots=Plots)
        s_avalancheanalysis.Avalanches(MovingAverage = False, T=N-N*0.4)
        s_avalancheanalysis.MomentAnalysis()
        s_avalancheanalysis.DMomentAnalysis()
        Moment[r] = np.array(s_avalancheanalysis.Moments_grads)
        DTau[r] = np.array([s_avalancheanalysis.D, s_avalancheanalysis.tau_s])
    Moment = np.rollaxis(Moment,1)
    DTau = np.rollaxis(DTau,1)
    Moment_data = np.array([s_avalancheanalysis.Moment_k,np.mean(Moment, 1), np.std(Moment, 1)/np.sqrt(R)])
    s_avalancheanalysis.DMomentAnalysisPlot(data=Moment_data)
    print "D: ", np.mean(DTau[0]), " +/- ", np.std(DTau[0])/np.sqrt(R)
    print "Tau_s: ", np.mean(DTau[1]), " +/- ", np.std(DTau[1])/np.sqrt(R)
     
def SavedHeightsRun(heights, Plots):
    """Produces all the graphs for the height analysis"""
    #s_heightanalysis = sha.SavedHeightsAnalysis(location = "", Plots=[[1000000, 8, 53], [1000000,16, 205], [1000000, 32, 865], [1000000,64, 3488], [1000000, 128, 14028], [1000000, 256, 56027], [1000000, 512, 225305]])
    s_heightanalysis = sha.SavedHeightsAnalysis(data = heights, Plots=Plots)
    s_heightanalysis.Heights(MovingAverage = False, CalculateCriticals = False, t0=100000, T=600000)
    s_heightanalysis.Gradients(points_excluded=0)
    s_heightanalysis.HeightsPlot()
    s_heightanalysis.LogHeightGradientPlot()
    s_heightanalysis.CrossOverTimePlot()
    s_heightanalysis.HeightOfPilePlot()
    s_heightanalysis.GradientOfPilePlot()
    s_heightanalysis.SDHeightOfPilePlot()
    s_heightanalysis.CorrectionToScalingHeightPlot(a0=1.8, optimise=True)
    s_heightanalysis.CorrectionToScalingSDHeightPlot(b0=1.8, alpha=3, optimise=True)
    s_heightanalysis.CollapsedHeightsPlot(D=None)
    s_heightanalysis.ProbabilitiesPlot()
    s_heightanalysis.CollapsedProbabilitiesPlot()
    
def SavedAvalanchesRun(avalanches, Plots):
    """Produces all the graphs for the avalanche-size analysis"""
    #s_avalancheanalysis = saa.SavedAvalancheAnalysis(location = "", Plots=[[1000000, 8, 53], [1000000,16, 205], [1000000, 32, 865], [1000000,64, 3488], [1000000, 128, 14028], [1000000, 256, 56027], [1000000, 512, 225305]])
    s_avalancheanalysis = saa.SavedAvalancheAnalysis(data = avalanches, Plots = Plots)
    s_avalancheanalysis.Avalanches(MovingAverage = False, T=700000)
    s_avalancheanalysis.DandTauEstimate(a=1.3)
    s_avalancheanalysis.RawLogBinAndTauS(a=1.4)
    s_avalancheanalysis.ProbabilitiesPlot(a=1.4)    
    s_avalancheanalysis.MomentAnalysis()
    s_avalancheanalysis.MomentAnalysisPlot()
    s_avalancheanalysis.DMomentAnalysisPlot()
    s_avalancheanalysis.CollapsedProbabilitiesPlot(a=1.5, D=None, tau_s=None)
    
def SavedDropsRun(drops, Plots):
    """Produces all the graphs for the drop size analysis"""
    #s_dropanalysis = sda.SavedDropAnalysis(location = "", Plots=[[1000000, 8, 53], [1000000,16, 205], [1000000, 32, 865], [1000000,64, 3488], [1000000, 128, 14028], [1000000, 256, 56027], [1000000, 512, 225305]])
    s_dropanalysis = sda.SavedDropAnalysis(data = drops, Plots = Plots)
    s_dropanalysis.Drops(MovingAverage = False)
    s_dropanalysis.RawLogBinAndTauS(a=1.6)
    s_dropanalysis.ProbabilitiesPlot(a=1.6)
    s_dropanalysis.CollapsedProbabilitiesLPlot(a=1.6, alpha=1.2, D=1.25)
    s_dropanalysis.CollapsedProbabilitiesDropPlot(a=1.6, tau_d=1.0, D=1.25)
    

#OsloModelPython(N=10000, L=8, p=0.5)
#Plots, heights, avalanches, drops = MutliOsloModelCPP(N=1000000, Lmin=3, Lmax=9, p=0.5)