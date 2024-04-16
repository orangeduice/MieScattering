# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 22:38:37 2024

@author: osjac
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import scipy
import matplotlib.gridspec as gridspec




def psi_e(phi, theta, m, n):
    return np.cos(m*phi)*scipy.special.lpmv(m,n,np.cos(theta))


def psi_o(phi, theta, m, n):
    return np.sin(m*phi)*scipy.special.lpmv(m,n,np.cos(theta))


def polar2cart(phi, theta, r):
    x = r * np.sin(phi) * np.cos(theta)
    y = r *np.sin(theta) * np.sin(phi)
    z = r *np.cos(phi)
    return x, y, z


def sphercialHarmCart(m,n):
    phi = np.linspace(0, 2*np.pi,num=100)
    theta = np.linspace(0, np.pi,num=100)
    
    PHI, THETA = np.meshgrid(phi, theta)
    xyz = np.array([np.sin(THETA) * np.sin(PHI),
                    np.sin(THETA) * np.cos(PHI),
                    np.cos(THETA)])
    Y = psi_e(PHI,THETA,m,n)
    Yx, Yy, Yz = np.abs(Y) * xyz

    return Yx, Yy, Yz, Y


def surfacePlotCart(ax, x, y, z, colorMap = 'coolwarm', colorScale = None):
    max_x = np.max(x)
    max_y = np.max(y)
    max_z = np.max(z)
    MAX = max([max_x,max_y,max_z]).round()
    ax.plot([-MAX*0.10, MAX*1.25], [0,0], [0,0], c='black', lw=1,zorder=10)
    ax.plot([0,0], [-MAX*0.10, MAX*1.25], [0,0], c='black', lw=1,zorder=10)
    ax.plot([0,0], [0,0], [-MAX*0.10, MAX*1.25], c='black', lw=1,zorder=10)

    if MAX != 0.0:
        ax.text(MAX*1.5, 0, 0, "x")
        ax.text(0, MAX*1.5, 0, "y")
        ax.text(0, 0, MAX*1.5, "z")

    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('coolwarm'))
    cmap.set_clim(-1, 1)
    
    ax.plot_surface(x, y, z, facecolors=cmap.to_rgba(colorScale), rstride=2, cstride=2)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-MAX,MAX)
    ax.set_ylim(-MAX,MAX)
    ax.set_zlim(-MAX,MAX)
    ax.set_box_aspect([1,1,1])
    ax.axis('off')
    #plt.show()


def multiSphercialHarmPlot(fig,m_max):
    spec = gridspec.GridSpec(ncols=m_max+1, nrows=m_max+1, figure=fig)
    for n in range(0,m_max+1):
        for m in range(0,m_max+1):
            print(m,n)
            ax = fig.add_subplot(spec[n, m], projection='3d')
            x, y, z, vaules = sphercialHarmCart(m,n)
            surfacePlotCart(ax, x, y, z, colorScale = vaules.real)




if __name__ == '__main__':
    fig = plt.figure(figsize=(10, 10), dpi=200)

    #Plot single
    #m, n = 0, 3
    #ax = fig.add_subplot(projection='3d')
    #x, y, z, vaules = sphercialHarmCart(m,n)
    #surfacePlotCart(x, y, z, colorScale = vaules.real)

    #Plot multiple
    m_max = 4
    multiSphercialHarmPlot(fig, m_max)

    plt.tight_layout()
    plt.savefig('sph_harm.png')
    plt.show()