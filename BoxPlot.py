import numpy as np
import matplotlib.pyplot as plt

class BoxPlot:
    def __init__(self, x_max, y_max):
        self.x_max = x_max
        self.y_max= y_max
        self.plot = 0
        
    def Start(self):
        self.plot = 1
        plt.figure('Animation')
        plt.title('Oslo Model Animation')
        plt.scatter(0, 0, s=20**2,marker='s', edgecolors='none')
        plt.axis([0, self.x_max, 0, self.y_max])
        plt.show()
            
    def Iterate(self, system):
        if self.plot == 0:
            return None
        system_plot = []
        L = len(system)-1
        for i in range(L):
            for y in range(int(system[i][0])):
                system_plot.append([i,y])
        system_plot = np.array(system_plot)
        plt.clf()
        plt.title('Oslo Model Animation')
        plt.axis([0, self.x_max, 0, self.y_max])
        plt.scatter(system_plot[:,0], system_plot[:,1], s=20**2,marker='s', edgecolors='none')
        plt.pause(0.0001)