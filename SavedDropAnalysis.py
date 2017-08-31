import numpy as np
import matplotlib.pyplot as plt
import DropAnalysis as ha
        
class SavedDropAnalysis:
    def __init__(self, location = "", data = [], Plots=[]):
        self.location = location
        self.data = data
        self.Plots = np.array(Plots)[:,:2]
        self.dropanalysis = []
        self.criticals = np.array(Plots)[:,2]
        self.D = self.tau_s = 0
        
    def Drops(self, MovingAverage = False):
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            
            if self.location != "":
                Drop = np.loadtxt(self.location + "Iter_" + str(N) + "_L_" + str(L))
            elif self.data != []:
                Drop = self.data[i]
            else:
                raise Exception("No location or data")
                
            self.dropanalysis.append(ha.DropAnalysis(N=N, L=L, Data=Drop))
            
            if MovingAverage: self.dropanalysis[i].MovingAverage(W=25)
            
            self.dropanalysis[i].FindCrossOverTime(self.criticals[i])
            self.dropanalysis[i].Probabilities()
            
            print "L:", str(L)
        
    def RawLogBinAndTauS(self, a):
        self.dropanalysis[-1].MovingAverage(W=0)
        self.dropanalysis[-1].Probabilities()
        self.dropanalysis[-1].LogBinProbabilities(a)
        self.dropanalysis[-1].ProbabilitiesVsLogBinPlot()
        
    def ProbabilitiesPlot(self, a):
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.dropanalysis[i].MovingAverage(W=0)
            self.dropanalysis[i].LogBinProbabilities(a)
            self.dropanalysis[i].LogBinProbabilitiesPlot()
            
    def CollapsedProbabilitiesLPlot(self, a, alpha=None, D=None):
        if alpha != None: self.alpha = alpha
        if D != None: self.D = D
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.dropanalysis[i].MovingAverage(W=0)
            self.dropanalysis[i].LogBinProbabilities(a)
            self.dropanalysis[i].CollapsedLogBinProbabilitiesLPlot(self.alpha, self.D)
            
    def CollapsedProbabilitiesDropPlot(self, a, tau_d=None, D=None):
        if tau_d != None: self.tau_d = tau_d
        if D != None: self.D = D
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.dropanalysis[i].MovingAverage(W=0)
            self.dropanalysis[i].LogBinProbabilities(a)
            self.dropanalysis[i].CollapsedLogBinProbabilitiesDropPlot(self.tau_d, self.D)