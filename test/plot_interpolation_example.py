import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_theme(style = 'dark')

if __name__ == '__main__':
    ipts = np.loadtxt('interpolation_points.dat')
    pts = np.loadtxt('interpolated.dat')

    fig = plt.figure(
        figsize = (10, 8)
    )

    plt.plot(ipts[:, 0], ipts[:, 1], label = 'Interpolation points', linestyle = '--', marker = 'o',)
    plt.plot(pts[:, 0], pts[:, 1], label = 'Interpolated airfoil')

    plt.xlim((-0.1, 1.1))
    plt.ylim((-0.6, 0.6))

    plt.xlabel('x/c')
    plt.ylabel('y/c')
    
    plt.grid()

    plt.legend()

    plt.show()
