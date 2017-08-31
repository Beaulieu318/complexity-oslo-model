import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def HeightCorrectionToScaling(func):
    a0_range = np.linspace(1.73, 1.76, 100)
    
    Z = np.zeros((len(a0_range)))
    for a0_num in range(len(a0_range)):
        Z[a0_num] = func.CorrectionToScalingHeight(a0_range[a0_num])
           
    plt.figure('Height CTS')
    plt.plot(a0_range, Z)
    plt.xlabel('a0')
    plt.ylabel('1-r^2')
    plt.show()
    
def SDCorrectionToScaling2D(func):
    b0 = 50
    alpha_range = np.linspace(-2, 5, 10000)
    
    Z = np.zeros((len(alpha_range)))
    for alpha_num in range(len(alpha_range)):
        Z[alpha_num] = func.CorrectionToScalingSDHeight([b0, alpha_range[alpha_num]])
            
    plt.figure('SD CTS 2D')
    plt.plot(alpha_range, Z)
    plt.xlabel('alpha')
    plt.ylabel('1-r^2')
    plt.show()

def SDCorrectionToScaling3D(func):
    b0_range = np.linspace(0, 100, 100)
    alpha_range = np.linspace(-5, 5, 100)
    
    X, Y = np.meshgrid(b0_range, alpha_range)
    
    Z = np.zeros((X.shape[0], X.shape[1]))
    for b0_num in range(len(b0_range)):
        for alpha_num in range(len(alpha_range)):
            Z[b0_num][alpha_num] = func.CorrectionToScalingSDHeight([X[b0_num][alpha_num], Y[b0_num][alpha_num]])
            
    fig = plt.figure("SD CTS 3D")
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z)
    ax.set_xlabel('b0')
    ax.set_ylabel('alpha')
    ax.set_zlabel('1-r^2')
    plt.show()