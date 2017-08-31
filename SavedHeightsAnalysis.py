import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, optimize
import HeightsAnalysis as ha
        
class SavedHeightsAnalysis:
    """Analysis on the heights of the Oslo model"""
    def __init__(self, location = "", data = [], Plots=[]):
        self.location = location
        self.data = data
        self.Plots = np.array(Plots)[:,:2]
        self.heightanalysis = []
        self.criticals = np.array(Plots)[:,2]
        
    def Heights(self, MovingAverage = False, CalculateCriticals = True, t0=100000, T=1000000):
        """Required to run to analyse the data.
        Inputs: 
            MovingAverage - smooths the data using windowing
            CalculateCriticals - not recommended on big data because it takes a while
            t0 - the buffer after the critical time
            T - the maximum number of drives after the ciritcal time
            """
        self.COTs = []
        self.HOPs = []
        self.LHGs = []
        self.SDHOPs = []
        self.Moments = []
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            
            if self.location != "":
                heights = np.loadtxt(self.location + "Iter_" + str(N) + "_L_" + str(L))
            elif self.data != []:
                heights = self.data[i]
            else:
                raise Exception("No location or data")
                
            self.heightanalysis.append(ha.HeightsAnalysis(N=N, L=L, Data=heights))
            if MovingAverage: self.heightanalysis[i].MovingAverage(W=25)
            
            critical = None
            if not CalculateCriticals: critical = self.criticals[i]
            COT = self.heightanalysis[i].FindCrossOverTime(critical)[0]
            HOP = self.heightanalysis[i].FindHeightOfPile(t0=t0, T=T)
            SDHOP = self.heightanalysis[i].FindSDHeightOfPile(t0=t0, T=T)
            LHG = self.heightanalysis[i].FindLogHeightsGradient()
            
            self.COTs.append(COT)
            self.HOPs.append(HOP)
            self.SDHOPs.append(SDHOP)
            self.LHGs.append(LHG)
            
            print "L:", str(L), ", CrossOverTime:", COT, ", Steady state height:", HOP, ' +/- ', SDHOP
            
    def Gradients(self, points_excluded=0):
        """Required for most of the following graphs
        Input:
            points_excluded which ignores the inital points to show corrections to scaling and gives more accurate values"""
        self.COTs_grad, self.COTs_intercept = stats.linregress(np.log10(self.Plots[:,1][-points_excluded:]),np.log10(self.COTs[-points_excluded:]))[:2]
        self.HOPs_grad, self.HOPs_intercept = stats.linregress(self.Plots[:,1][-points_excluded:],self.HOPs[-points_excluded:])[:2]
        self.SDHOPs_grad, self.SDHOPs_intercept = stats.linregress(np.log10(self.Plots[:,1][-points_excluded:]),np.log10(self.SDHOPs[-points_excluded:]))[:2]
        self.LHG_grad, self.LHG_intercept = self.LHGs[-1]
        self.AvG = self.HOPs_grad
        print "Height Of Pile Gradient: ", self.HOPs_grad
        print "Uncertainty Of Height Of Pile Exponent of L: ", self.SDHOPs_grad
        print "Cross Over Time Exponent of L (D): ", self.COTs_grad
        print "Exponent for the transient configuration: ", self.LHG_grad
        
    def LogHeightGradientPlot(self, plot="log"):
        """Used in conjunction with HeightsPlot to show the line of best fit in the recurrent configuration
        Inputs: 
            the type of graph (log or linear)"""
        plt.figure('Height')
        x_range = np.logspace(0, 7, 100000)
        y_range = [10**self.LHG_intercept*x**self.LHG_grad for x in x_range]
        plt.figure('Height')
        plt.plot(x_range, y_range, 'k--', label='Line of best fit')
        if plot == "linear":
            plt.legend()
        elif plot == "log":
            plt.legend(loc=2)
        plt.show()
            
    def CrossOverTimePlot(self):
        """Plots the cross-over time with respect to system size"""
        x_range = np.linspace(0, self.Plots[:,1][-1])
        y_range = [10**self.COTs_intercept*x**self.COTs_grad for x in x_range]
        plt.figure('CrossOverTime')
        plt.scatter(self.Plots[:,1], self.COTs)
        plt.plot(x_range, y_range, 'k--', label='Line of best fit')
        plt.title('The cross-over time for the to reach stead state \n against the system size')
        plt.ylabel('Cross-over Time (t)')
        plt.xlabel('System Size')
        plt.grid(True)
        plt.legend()
        plt.show()
        
    def HeightOfPilePlot(self):
        x_range = np.linspace(0, self.Plots[:,1][-1])
        y_range = [self.HOPs_grad*x + self.HOPs_intercept for x in x_range]
        plt.figure('HeightOfPile')
        plt.errorbar(self.Plots[:,1], self.HOPs, yerr=self.SDHOPs, fmt='o')
        plt.plot(x_range, y_range, 'k--', label='Line of best fit')
        plt.title('The average height of the pile at stead state \n against the system size')
        plt.ylabel('Average height (<h(t; L)>)')
        plt.xlabel('System Size (L)')
        plt.grid(True)
        plt.legend()
        plt.show()
        
    def GradientOfPilePlot(self):
        """Plots the gradient of the pile with respect to system size (also the correction to scaling in the height"""
        plt.figure('GradientOfPile')
        x_range = np.linspace(0, self.Plots[:,1][-1])
        y_range = [self.AvG for x in x_range]
        plt.plot(x_range, y_range, 'k--', label='Expected average gradient')
        plt.legend()
        plt.scatter(self.Plots[:,1], self.HOPs/self.Plots[:,1])
        plt.title('The average gradient of the pile at stead state \n against the system size')
        plt.ylabel(r'<h(t; L)>/L')
        plt.xlabel('System Size (L)')
        plt.grid(True)
        plt.show()
        
    def SDHeightOfPilePlot(self, alpha=0, data=[]):
        """Plots the standard deviation of the height of the pile with respect to system size
        Inputs:
            can vary alpha the exponent of L
            can import data and show the uncertainties in the standard deviaiton."""
        plt.figure('SDOfHeightOfPile')
        if data == []:
            x_vals = self.Plots[:,1]
            y_vals = self.SDHOPs/self.Plots[:,1]**alpha
            yerr = None
            x_range = [0, self.Plots[:,1][-1]]
            y_range = [10**self.SDHOPs_intercept, 10**self.SDHOPs_intercept]
            plt.plot(x_range, y_range, 'r--', label='Expected standard deviation')
        else:
            x_vals = data[0]
            y_vals = data[1]/self.Plots[:,1]**alpha
            yerr = data[2]/self.Plots[:,1]**alpha
            self.SDHOPs_grad, self.SDHOPs_intercept = stats.linregress(np.log10(x_vals[-4:]),np.log10(y_vals[-4:]))[:2]
            x_range = [0, self.Plots[:,1][-1]]
            y_range = [10**self.SDHOPs_intercept, 10**self.SDHOPs_intercept]
            plt.plot(x_range, y_range, 'r--', label='Expected standard deviation')
        if alpha == 0:
            x_range = np.linspace(1, self.Plots[:,1][-1])
            y_range = [10**self.SDHOPs_intercept*x**self.SDHOPs_grad for x in x_range]
            plt.plot(x_range, y_range, 'k--', label='Line of best fit')
        plt.errorbar(x_vals, y_vals, yerr=yerr, fmt='o', label = "alpha=" + str(alpha))
        plt.title('The uncertainty of the average height of the pile at stead state \n against the system size')
        plt.ylabel(r'$\sigma/{L^\alpha}$')
        plt.xlabel('System Size (L)')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    def CorrectionToScalingHeight(self, a0):
        """Function to output the 1-R-squared value of the correction to scaling in the height"""
        ctsh_y = 1-self.HOPs/(a0*self.Plots[:,1])
        if np.any(ctsh_y < 0) : return 1
        self.CTSHs_grad, self.CTSHs_intercept, rval, pval, stderr = stats.linregress(np.log10(self.Plots[:,1]),np.log10(ctsh_y))
        return 1 - rval**2
        
    def CorrectionToScalingHeightPlot(self, a0, optimise=True):
        """Find a0 and omega1 in the correction to scaling of the heights
        Inputs:
            a0 is a trival estimation"""
        if optimise: a0 = optimize.minimize(self.CorrectionToScalingHeight, a0)['x'][0]
        self.CorrectionToScalingHeight(a0) 
        ctsh_y = (1-self.HOPs/(a0*self.Plots[:,1]))
        x_range = np.linspace(0, np.log10(self.Plots[:,1][-1]), 1000)
        y_range = [self.CTSHs_intercept + self.CTSHs_grad*x for x in x_range]
        plt.figure('CorrectionToScalingHeight')
        plt.scatter(np.log10(self.Plots[:,1]), np.log10(ctsh_y), label='Raw Data:\n a0=' + str(a0) + "\n omega1=" + str(-self.CTSHs_grad))
        plt.plot(x_range, y_range, 'k--', label='Line of best fit')
        plt.title('Loglog graph showing the correction to scaling \n of the height')
        plt.ylabel('1 - <h(t; L)>/(a0 x L)')
        plt.xlabel('System Size (L)')
        plt.grid(True)
        plt.legend()
        plt.show()
        self.a0 = a0
        self.AvG = a0
        print "Correction to Scaling Height - a0:", self.AvG
        print "Correction to Scaling Height - w1: ", -self.CTSHs_grad
        
    def CorrectionToScalingSDHeight(self, variables):
        b0 = variables[0]
        alpha = variables[1]
        ctsh_y = 1-self.SDHOPs/(b0*self.Plots[:,1]**alpha)
        if np.any(ctsh_y < 0) : return 1
        self.CTSSDHs_grad, self.CTSSDHs_intercept, rval, pval, stderr = stats.linregress(np.log10(self.Plots[:,1]),np.log10(ctsh_y))
        return 1 - rval**2
        
    def CorrectionToScalingSDHeightPlot(self, b0, alpha, optimise=True):
        if optimise: b0, alpha = optimize.minimize(self.CorrectionToScalingSDHeight, [b0, alpha])['x']
        self.CorrectionToScalingSDHeight([b0, alpha]) 
        ctsh_y = (1-self.SDHOPs/(b0*self.Plots[:,1]**alpha))
        x_range = np.linspace(0, np.log10(self.Plots[:,1][-1]), 1000)
        y_range = [self.CTSSDHs_intercept + self.CTSSDHs_grad*x for x in x_range]
        plt.figure('CorrectionToScalingSDHeight')
        plt.scatter(np.log10(self.Plots[:,1]), np.log10(ctsh_y), label='Raw Data:\n b0=' + str(b0) + "\n alpha=" + str(alpha) + "\n omega1=" + str(-self.CTSSDHs_grad))
        plt.plot(x_range, y_range, 'k--', label='Line of best fit')
        plt.title('Loglog graph showing the correction to scaling \n of the uncertainty of the height')
        plt.ylabel(r'$1 - <\sigma(t; L)>/(b_0 * L^\alpha)$')
        plt.xlabel('System Size (L)')
        plt.grid(True)
        plt.legend()
        plt.show()
        print "Correction to Scaling SD Height - b0: ", b0
        print "Correction to Scaling SD Height - alpha: ", alpha
        print "Correction to Scaling SD Height - gamma1: ", -self.CTSSDHs_grad
        
    def HeightsPlot(self, plot="log"):
        """Plots the heights against time"""
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.heightanalysis[i].MovingAverage(W=0)
            self.heightanalysis[i].Plot(plot=plot)
        
    def CollapsedHeightsPlot(self, D=None):
        """Plots the collapsed height against time
        Inputs:
            D can be varied"""
        if D != None: self.COTs_grad = np.array([D])
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.heightanalysis[i].MovingAverage(W=25)
            self.heightanalysis[i].CollapsedPlot(D=self.COTs_grad)
            
    def ProbabilitiesPlot(self):
        """Plots the probabilities against the height in the recurrent region."""
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.heightanalysis[i].MovingAverage(W=0)
            self.heightanalysis[i].Probabilities()
            self.heightanalysis[i].ProbabilitiesPlot()
            
    def CollapsedProbabilitiesPlot(self, plot=''):
        """Plots the collapsed proababilities against the height in the recurrent region.
        Inputs:
            plot = functional -> functional form of the standard deviation and the average height instead of the exact values"""
        for i in range(len(self.Plots)):
            N, L = self.Plots[i]
            self.heightanalysis[i].MovingAverage(W=0)
            self.heightanalysis[i].Probabilities()
            if plot == 'functional':
                av_h = self.AvG*self.Plots[:,1][i]
                st_dev = 10**self.SDHOPs_intercept*self.Plots[:,1][i]**self.SDHOPs_grad
            else:
                av_h=self.HOPs[i]
                st_dev=self.SDHOPs[i]
            self.heightanalysis[i].CollapsedProbabilitiesPlot(av_h, st_dev)