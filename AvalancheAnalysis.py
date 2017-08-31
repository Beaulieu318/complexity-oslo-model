import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import log_bin_CN_2016 as lb

class AvalancheAnalysis:
    def __init__(self, N=0, L=0, T=700000, Data=[]):
        self.N = N
        self.L = L
        self.T = T
        self.Data_original = np.array(Data)
        self.Data = np.array(Data)
        self.W = 0
        self.t_range_original = np.array(range(self.N))
        self.t_range = np.array(range(self.N))
        self.cot_x = self.cot_y = 0
        self.LBLOBF_xrange = self.LBLOBF_yrange = np.array([])
        
    def SaveData(self, location=""):
        np.savetxt(location+"/Iter_"+str(self.N)+"_L_"+str(self.L), self.Data)
        
    def MovingAverage(self, W):
        if W == 0:
            self.Data = self.Data_original
            self.t_range = self.t_range_original
        else:
            self.W = W
            avalanche_MA = []
            for t in range(len(self.Data)):
                avalanche_MA.append(1.0/(2.0 * float(W) + 1.0) * np.sum(self.Data[t-W:t+W]))
            self.Data = np.array(avalanche_MA)[W:-W]
            self.t_range = self.t_range[W:-W]
        
    def FindCrossOverTime(self, cot=0):
        self.cot_x = cot
        self.cot_y = self.Data[cot]
        return self.cot_x, self.cot_y
                
    def Probabilities(self):
        data = self.Data[self.cot_x:self.cot_x+self.T].astype(int)
        self.avalanche = np.arange(min(data), max(data)+1)
        count = np.bincount(data)[min(data): max(data)+1]
        self.prob = count.astype(float)/np.sum(count)
        self.avalanche = self.avalanche[1:]
        self.prob = self.prob[1:]
        
    def LogBinProbabilities(self, a=1.4):
        self.avalanche_bin, self.prob_bin = lb.log_bin(self.Data[self.cot_x:self.cot_x+self.T], bin_start=0., first_bin_width=1., a=a, datatype='float', drop_zeros=False, debug_mode=False)
        self.avalanche_bin = np.array(self.avalanche_bin[1:])
        self.prob_bin = np.array(self.prob_bin[1:])
        
    def Moment(self, k):
        data = self.Data[self.cot_x:self.cot_x+self.T]
        return np.mean(np.float64(data)**float(k))
        
    def LogBinLineOfBestFit(self):
        self.LBLOBF_grad, self.LBLOBF_intercept = stats.linregress(np.log10(self.avalanche_bin[:10000]),np.log10(self.prob_bin[:10000]))[:2]
        self.LBLOBF_xrange = np.logspace(0, 6, 50)
        self.LBLOBF_yrange = np.array([10**self.LBLOBF_intercept*x**self.LBLOBF_grad for x in self.LBLOBF_xrange])
        return self.LBLOBF_grad
        
    def Plot(self):
        avalanche = self.Data
        plt.figure('AvalancheSize')
        plt.title('Avalanche-Size')
        plt.plot(self.t_range, avalanche, label = 'L=' + str(self.L))
        plt.plot((self.cot_x, self.cot_x), (0, self.cot_y), 'k-')
        plt.plot((0, self.cot_x), (self.cot_y, self.cot_y), 'k-')
        plt.ylabel('s(t; L) (Avalanche-size)')
        plt.xlabel('t (Number of blocks added)')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()
        
    def ProbabilitiesVsLogBinPlot(self):
        plt.figure('ProbabilitiesVsLogBinPlot')
        plt.title('Avalanche-size probabilities')
        plt.scatter(np.log10(self.avalanche), np.log10(self.prob), label = 'L=' + str(self.L) + ": raw data")
        plt.plot(np.log10(self.avalanche_bin), np.log10(self.prob_bin), label = 'L=' + str(self.L) + ": logbin")
        plt.plot(np.log10(self.LBLOBF_xrange), np.log10(self.LBLOBF_yrange), 'r--', label = "Line of best fit")
        plt.ylabel('log(P(s; L))')
        plt.xlabel('log(s)')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    def LogBinProbabilitiesPlot(self):
        plt.figure('AvalancheProbabilitiesPlot')
        plt.title('Avalanche-size probabilities')
        plt.plot(self.avalanche_bin, self.prob_bin, label = 'L=' + str(self.L))
        plt.ylabel('P(s; L)')
        plt.xlabel('s')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    def CollapsedLogBinProbabilitiesPlot(self, tau_s, D):
        plt.figure('CollapsedAvalancheProbabilitiesPlot')
        plt.title('Collapsed avalanche-size probabilities')
        plt.plot(self.avalanche_bin/np.array([self.L])**float(D), self.prob_bin*self.avalanche_bin**float(tau_s), label = 'L=' + str(self.L))
        plt.ylabel(r'$s^{\tau_s}$P(s; L)')
        plt.xlabel(r's/$L^D$')
        plt.legend(loc=3)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()