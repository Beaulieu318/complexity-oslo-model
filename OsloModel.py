import numpy as np
import BoxPlot as bp
import HeightsAnalysis as ha

class OsloModel:
    def __init__(self,L):
        self.L = L
        self.system = np.zeros([self.L+1,2])
        self.plot = 0
        self.heights = []
        self.critical = 0
        self.foundcritical = False
            
    def Drive(self):
        self.system[0][0] += 1
        
    def Relaxation(self, n, p):
        self.Box_Plot.Iterate(self.system)
        relaxed = False
        s = 0
        d = 0
        while not relaxed:
            relaxed = True
            for i in range(self.L):
                z = self.system[i][0] - self.system[i+1][0]
                if z > self.system[i][1]:
                    if i < self.L -1:
                        self.system[i][0] -= 1
                        self.system[i+1][0] += 1
                    else:
                        self.system[i][0] -= 1
                        self.CriticalN(n)
                        d += 1
                    self.system[i][1] = np.random.choice(np.arange(1, 3), p=[p, 1-p])
                    self.Box_Plot.Iterate(self.system)
                    s += 1
                    relaxed = False
        self.avalanchesize[n] = s
        self.dropsize[n] = d
                
    def Iteration(self, N, p, Box=True):
        
        self.Box_Plot = bp.BoxPlot(x_max = 20, y_max = 40)
        if Box == True:
            self.Box_Plot.Start()
        
        self.avalanchesize = np.zeros([N, 1])
        self.dropsize = np.zeros([N, 1])
        
        for i in range(self.L):
            self.system[i][1] = np.random.choice(np.arange(1, 3), p=[p, 1-p])
        
        for n in range(N):
            self.Drive()
            self.Relaxation(n, p)
            self.Height(self.system)
            
        self.heights = np.array(self.heights)
        self.avalanchesize = self.avalanchesize.flatten()
        self.dropsize = self.dropsize.flatten()
            
    def CriticalN(self, n):
        if self.foundcritical == False:
            self.critical = n
            self.foundcritical = True
    
    def Height(self, system):
        self.heights.append(system[0][0])
        return np.array(self.heights)
        
    def AverageSlope(self, system):
        return self.Height/self.L