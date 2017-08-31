import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import AvalancheAnalysis as ha
        
class SavedAvalancheAnalysis:
    """Analysis on the avalanche-size probabilities of the Oslo model"""
    def __init__(self, location = "", data=[], Plots=[]):
        self.location = location
        self.data = data       
        self.Plots = np.array(Plots)[:,:2]
        self.avalancheanalysis = []
        self.criticals = np.array(Plots)[:,2]
        self.D = self.tau_s = 0
        
    def Avalanches(self, MovingAverage = False, T=700000):
        """Required to run analyse on the data
        Inputs:
            MovingAverage
            The number of avalanches included in the analysis."""
        self.Moment_k = range(1, 7)
        self.Moments = []
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            
            if self.location != "":
                Avalanche = np.loadtxt(self.location + "Iter_" + str(N) + "_L_" + str(L))
            elif self.data != []:
                Avalanche = self.data[i]
            else:
                raise Exception("No location or data")
                
            self.avalancheanalysis.append(ha.AvalancheAnalysis(N=N, L=L, T=T, Data=Avalanche))
            
            if MovingAverage: self.avalancheanalysis[i].MovingAverage(W=25)
            
            self.avalancheanalysis[i].FindCrossOverTime(self.criticals[i])
            self.avalancheanalysis[i].Probabilities()
            self.Moments.append([self.avalancheanalysis[i].Moment(j) for j in self.Moment_k])
                        
            print "L:", str(L)
        self.Moments = np.array(self.Moments)
        self.Moments = np.rollaxis(self.Moments, 1)
        
    def DandTauEstimate(self, a, points_excluded=0):
        """Finds a primitive estimate of D and Tau"""
        self.LastPoint = []
        for i in range(len(self.Plots)):
            self.avalancheanalysis[i].LogBinProbabilities(a)
            self.LastPoint.append(self.avalancheanalysis[i].avalanche_bin[-1])
        self.LPs_grad, self.LPs_intercept = stats.linregress(np.log10(self.Plots[:,1][-points_excluded:]),np.log10(self.LastPoint[-points_excluded:]))[:2]
        self.tau_s = -self.avalancheanalysis[-1].LogBinLineOfBestFit()
        self.D = self.LPs_grad
        print "Tau_s: ", self.tau_s
        print "D: ", self.D
        
    def MomentAnalysis(self):
        """Performs moment anlysis on the data"""
        self.Moments_grads = []
        self.Moments_intercepts = []
        for Moment in self.Moments:
            self.Moments_grad, self.Moments_intercept = stats.linregress(np.log10(self.Plots[:,1][-4:]),np.log10(Moment[-4:]))[:2]
            self.Moments_grads.append(self.Moments_grad)
            self.Moments_intercepts.append(self.Moments_intercept)
        
    def MomentAnalysisPlot(self):
        """Plots the moment analysis graph"""
        self.MomentAnalysis()
        k = 1
        for i in range(len(self.Moments)):
            x_range = np.linspace(1, self.Plots[:,1][-1], 10000)
            y_range = [10**self.Moments_intercepts[i]*x**self.Moments_grads[i] for x in x_range]
            plt.figure('MomentsAnalysis')
            plt.scatter(self.Plots[:,1], self.Moments[i])
            plt.plot(x_range, y_range, '--', label='k='+str(k))
            plt.title('kth Moment against the system size \n for the average avalanche-size')
            plt.ylabel(r'<$s^k$>')
            plt.xlabel('System Size (L)')
            plt.xscale('log')
            plt.yscale('log')
            plt.grid(True)
            plt.legend(loc=2)
            plt.show()
            k += 1
            
    def DMomentAnalysis(self):
        """Finds D and tau_s from moment analysis"""
        x_vals = self.Moment_k
        y_vals = self.Moments_grads
        k_grad, k_intercept = stats.linregress(x_vals, y_vals)[:2]
        self.D = k_grad
        self.tau_s = 1 - k_intercept/self.D
            
    def DMomentAnalysisPlot(self, data = []):
        """Plots the final graph from moment analysis"""
        self.MomentAnalysis()
        if data == []:
            x_vals = self.Moment_k
            y_vals = self.Moments_grads
            yerr = None
        else:
            x_vals = data[0]
            y_vals = data[1]
            yerr = data[2]
        plt.figure('DMomentsAnalysis')
        k_grad, k_intercept = stats.linregress(x_vals, y_vals)[:2]
        x_range = np.linspace(0, max(x_vals), 10000)
        y_range = [k_grad*x + k_intercept for x in x_range]
        plt.errorbar(x_vals, y_vals, yerr=yerr, fmt='o')
        plt.plot(x_range, y_range, 'k--')
        plt.title('Estimating D against the moment k \n for the average avalanche-size')
        plt.ylabel(r'D(1+k+$\tau_s$)')
        plt.xlabel('k')
        plt.grid(True)
        plt.show()
        self.D = k_grad
        self.tau_s = 1 - k_intercept/self.D
        print 'Moment Analysis - D: ', self.D
        print 'Moment Analysis - tau_s', self.tau_s
        
    def RawLogBinAndTauS(self, a):
        """Plots the raw data against the log binned data for the largest system size
        Inputs:
            the 'a' value used to increase the log bins"""
        self.avalancheanalysis[-1].MovingAverage(W=0)
        self.avalancheanalysis[-1].Probabilities()
        self.avalancheanalysis[-1].LogBinProbabilities(a)
        self.tau_s = -self.avalancheanalysis[-1].LogBinLineOfBestFit()
        self.avalancheanalysis[-1].ProbabilitiesVsLogBinPlot()
        print "Negative gradient of LogBin (tau_s) :", self.tau_s
        
    def ProbabilitiesPlot(self, a):
        """Plots the probabilities of all the system sizes on one graph
        Inputs:
            the 'a' value can be varied."""
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.avalancheanalysis[i].MovingAverage(W=0)
            self.avalancheanalysis[i].LogBinProbabilities(a)
            self.avalancheanalysis[i].LogBinProbabilitiesPlot()
            
    def CollapsedProbabilitiesPlot(self, a, tau_s=None, D=None):
        """Plots the collapsed probabilities of all the system sizes on one graph
        Inputs:
            the 'a' value can be varied.
            an approximation of tau_s and D can be given for testing different values."""
        if tau_s != None: self.tau_s = tau_s
        if D != None: self.D = D
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.avalancheanalysis[i].MovingAverage(W=0)
            self.avalancheanalysis[i].LogBinProbabilities(a)
            self.avalancheanalysis[i].CollapsedLogBinProbabilitiesPlot(self.tau_s, self.D)