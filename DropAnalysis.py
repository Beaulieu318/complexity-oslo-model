import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import log_bin_CN_2016 as lb

class DropAnalysis:
    def __init__(self, N=0, L=0, Data=[]):
        self.N = N
        self.L = L
        self.Data_original = np.array(Data)
        self.Data = np.array(Data)
        self.W = 0
        self.t_range_original = np.array(range(self.N))
        self.t_range = np.array(range(self.N))
        self.cot_x = self.cot_y = 0
        self.LBLOBF_xrange = self.LBLOBF_yrange = np.array([])
        self.T = 700000
        
    def SaveData(self, location = ""):
        np.savetxt(location + "/Iter_"+str(self.N)+"_L_"+str(self.L), self.Data)
        
    def MovingAverage(self, W):
        if W == 0:
            self.Data = self.Data_original
            self.t_range = self.t_range_original
        else:
            self.W = W
            drop_MA = []
            for t in range(len(self.Data)):
                drop_MA.append(1.0/(2.0 * float(W) + 1.0) * np.sum(self.Data[t-W:t+W]))
            self.Data = np.array(drop_MA)[W:-W]
            self.t_range = self.t_range[W:-W]
        
    def FindCrossOverTime(self, cot=0):
        self.cot_x = cot
        self.cot_y = self.Data[cot]
        return self.cot_x, self.cot_y
                
    def Probabilities(self):
        data = self.Data[self.cot_x:self.cot_x+self.T].astype(int)
        self.drop = np.arange(min(data), max(data)+1)[1:]
        count = np.bincount(data)[min(data): max(data)+1][1:]
        self.prob = count.astype(float)/np.sum(count)
        self.drop = self.drop
        self.prob = self.prob
        
    def LogBinProbabilities(self, a=1.4):
        self.drop_bin, self.prob_bin = lb.log_bin(self.Data[self.cot_x:self.cot_x+self.T], bin_start=0., first_bin_width=1., a=a, datatype='float', drop_zeros=True, debug_mode=False)
        self.drop_bin = np.array(self.drop_bin[1:])
        self.prob_bin = np.array(self.prob_bin[1:])
        
    def Moment(self, k):
        data = self.Data[self.cot_x:self.cot_x+self.T]
        return np.mean(np.float64(data)**float(k))
        
    def LogBinLineOfBestFit(self):
        self.LBLOBF_grad, self.LBLOBF_intercept = stats.linregress(np.log10(self.drop_bin[:10000]),np.log10(self.prob_bin[:10000]))[:2]
        self.LBLOBF_xrange = np.logspace(0, 6, 50)
        self.LBLOBF_yrange = np.array([10**self.LBLOBF_intercept*x**self.LBLOBF_grad for x in self.LBLOBF_xrange])
        return self.LBLOBF_grad
        
    def Plot(self):
        drop = self.Data
        plt.figure('DropSize')
        plt.title('Oslo model Drop-Size')
        plt.plot(self.t_range, drop, label = 'L=' + str(self.L))
        plt.plot((self.cot_x, self.cot_x), (0, self.cot_y), 'k-')
        plt.plot((0, self.cot_x), (self.cot_y, self.cot_y), 'k-')
        plt.ylabel('d(t; L) (Drop-size)')
        plt.xlabel('t (Number of blocks added)')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()
        
    def ProbabilitiesVsLogBinPlot(self):
        plt.figure('ProbabilitiesVsLogBinPlot')
        plt.title('Drop probability of Oslo model')
        plt.scatter(np.log10(self.drop), np.log10(self.prob), label = 'L=' + str(self.L) + ": raw data")
        plt.plot(np.log10(self.drop_bin), np.log10(self.prob_bin), label = 'L=' + str(self.L) + ": logbin")
        plt.plot(np.log10(self.LBLOBF_xrange), np.log10(self.LBLOBF_yrange), 'r--', label = "Line of best fit")
        plt.ylabel('P(d; L)')
        plt.xlabel('d')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    def LogBinProbabilitiesPlot(self):
        plt.figure('DropProbabilitiesPlot')
        plt.title('Drop probability of Oslo model')
        plt.plot(self.drop_bin, self.prob_bin, label = 'L=' + str(self.L))
        plt.ylabel('P(d; L)')
        plt.xlabel('d')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    def CollapsedLogBinProbabilitiesLPlot(self, alpha, D):
        plt.figure('CollapsedDropProbabilitiesLPlot')
        plt.title('Drop probability of Oslo model')
        plt.plot(self.drop_bin/np.array([self.L])**float(D), self.prob_bin*np.array([self.L])**float(alpha), label = 'L=' + str(self.L))
        plt.ylabel(r'$L^{\alpha}$P(d; L)')
        plt.xlabel(r'd/$L^D$')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()
        
    def CollapsedLogBinProbabilitiesDropPlot(self, tau_d, D):
        plt.figure('CollapsedDropProbabilitiesDropPlot')
        plt.title('Drop probability of Oslo model')
        plt.plot(self.drop_bin/np.array([self.L])**float(D), self.prob_bin*self.drop_bin**float(tau_d), label = 'L=' + str(self.L))
        plt.ylabel(r'$d^{\tau_d}$P(d; L)')
        plt.xlabel(r'd/$L^D$')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()