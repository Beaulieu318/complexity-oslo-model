import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

class HeightsAnalysis:
    def __init__(self, N=0, L=0, Data=[]):
        self.N = N
        self.L = L
        self.Data_original = np.array(Data)
        self.Data = np.array(Data)
        self.W = 0
        self.t_range_original = np.array(range(self.N))
        self.t_range = np.array(range(self.N))
        self.cot_x = self.cot_y = 0
        
    def SaveData(self, location=""):
        np.savetxt(location + "/Iter_"+str(self.N)+"_L_"+str(self.L), self.Data)
        
    def MovingAverage(self, W=0):
        if W == 0:
            self.Data = self.Data_original
            self.t_range = self.t_range_original
        else:
            self.W = W
            heights_MA = []
            for t in range(len(self.Data)):
                heights_MA.append(1.0/(2.0 * float(W) + 1.0) * np.sum(self.Data[t-W:t+W]))
            self.Data = np.array(heights_MA)[W:-W]
            self.t_range = self.t_range[W:-W]
        
    def FindCrossOverTime(self, cot=None):
        if cot != None:
            self.cot_x = cot
            self.cot_y = self.Data[cot]
            return self.cot_x, self.cot_y
        V = 25*self.L
        for i in range(1, len(self.Data)-V):
            gradients = stats.linregress(np.log10(self.t_range[i: i+V]),np.log10(self.Data[i: i+V]))[0]
            if gradients < 0.0001:
                self.cot_x = self.W+i
                self.cot_y = self.Data[i]
                return self.cot_x, self.cot_y
                
    def FindHeightOfPile(self, t0=100000, T=1000000):
        return np.mean(self.Data[self.cot_x+t0:self.cot_x+t0+T])
        
    def FindSDHeightOfPile(self, t0=100000, T=100000):
        return np.std(self.Data[self.cot_x+t0:self.cot_x+t0+T])
        
    def FindLogHeightsGradient(self):
        return stats.linregress(np.log10(self.t_range[1: self.cot_x]),np.log10(self.Data[1: self.cot_x]))[:2]
        
    def Probabilities(self):
        data = self.Data[self.cot_x:].astype(int)
        self.height_prob = np.arange(min(data), max(data)+1)
        count = np.bincount(data)[min(data): max(data)+1]
        self.prob = count.astype(float)/np.sum(count)
        
    def Plot(self, plot="log"):
        heights = self.Data
        plt.figure('Height')
        plt.title('The heights')
        plt.plot(self.t_range, heights, label = 'L=' + str(self.L))
        plt.plot((self.cot_x, self.cot_x), (0, self.cot_y), 'k-')
        plt.plot((0, self.cot_x), (self.cot_y, self.cot_y), 'k-')
        plt.ylabel('h(t; L) (Number of blocks at i=0)')
        plt.xlabel('t (Number of blocks added)')
        if plot == "linear":
            plt.legend()
        elif plot == "log":
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=2)
        plt.grid(True)
        plt.show()
        
    def CollapsedPlot(self, D):
        heights = self.Data
        plt.figure('Collapsed Height')
        plt.title('The collapsed heights')
        plt.plot(self.t_range/float(self.L)**D, heights/float(self.L), label = 'L=' + str(self.L))
        plt.ylabel(r'h(t; L) / L')
        plt.xlabel(r'$t / L^D$')
        plt.legend(loc=2)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()
        
    def ProbabilitiesPlot(self):
        plt.figure('Probabilities')
        plt.title('The height probabilities')
        plt.plot(self.height_prob, self.prob, label = 'L=' + str(self.L))
        plt.ylabel('P(h; L)')
        plt.xlabel('h')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    def CollapsedProbabilitiesPlot(self, av_h, st_dev):
        plt.figure('Collapsed Probabilities')
        plt.title('The collapsed height probabilities')
        plt.plot((self.height_prob-av_h)/st_dev, self.prob*st_dev, label = 'L=' + str(self.L))
        plt.ylabel(r'$P(h; L)/\sigma$')
        plt.xlabel(r'$(h-<h>)/\sigma$')
        plt.legend()
        plt.grid(True)
        plt.show()