# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 17:46:27 2024

@author: okcj1g19
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
plt.style.use('default')

def Riccati_1(n, x):
    jn = special.spherical_jn(n, x)
    jn_dir = special.spherical_jn(n, x, derivative=True)
    r0 = x*jn
    r1 = jn + x*jn_dir
    return r0, r1


def Riccati_2(n, x):
    jn = special.spherical_jn(n, x)
    jn_dir = special.spherical_jn(n, x, derivative=True)
    yn = special.spherical_yn(n, x)
    yn_dir = special.spherical_yn(n, x, derivative=True)
    hn = jn + 1j*yn
    hnp = jn_dir + 1j*yn_dir
    r0 = x*hn
    r1 = hn + x*hnp
    return r0, r1


def scatCoff(n, m, x):
    if m == 1:
        return x*0, x*0
    jn, jn_p = Riccati_1(n, x)
    yn, yn_p = Riccati_2(n, x)
    jnm, jnm_p = Riccati_1(n, m*x)
    a = (m*jnm*jn_p - jn*jnm_p)/(m*jnm*yn_p - yn*jnm_p)
    b = (jnm*jn_p - m*jn*jnm_p)/(jnm*yn_p - m*yn*jnm_p)
    return a, b


def ScatteringCrossSection(m, a, i_var, f_var, res = 1000, norm = True, n_max = 10, multipole = False, prefix = "nm", sizeParameter = False):
    prefixes =  {"nm":10**-9, "um":10**-6, "mm":10**-3, "m":1}
    pfix = prefixes[prefix]

    a = a*pfix


    if not sizeParameter:
        i_var = i_var*pfix
        f_var = f_var*pfix
        i_x = (2*np.pi*a)/(i_var)
        f_x = (2*np.pi*a)/(f_var)
        x = np.linspace(i_x, f_x, num = res)
        k = x/a
        output_var = ((2*np.pi)/k)
    else:
        x = np.linspace(i_var, f_var, num = res)
        k = x/a
        output_var = x


    #x = np.linspace(x_i, x_f,num = res)
    #x = np.logspace(x_i, x_f,num = res,base=2)
    #x = np.geomspace(x_i, x_f,num = res)



    consts = (2*np.pi)/np.square(k)



    if norm:
        scaleFator = a**2
    else:
        scaleFator = 1
    summy = 0
    multipoles = []
    for n in range(1,n_max+1):
        aCoeff, bCoeff = scatCoff(n, m, x)
        aCoeff_norm = np.square(np.abs(aCoeff))
        bCoeff_norm = np.square(np.abs(bCoeff))
        summy = summy + (2*n + 1) * (aCoeff_norm + bCoeff_norm)
        if multipole:
            multipoles.append([(consts *(2*n + 1)*aCoeff_norm)/scaleFator, (consts *(2*n + 1)*bCoeff_norm)/scaleFator])
    scattering = (consts * summy)/scaleFator
    return output_var, scattering, multipoles


def PlotScatteringCrossSection(ax, waveLengths, scattering, max_dis = 2, multipoles = None, norm = True, sizeParameter = False, prefix = "nm"):
    keys = [["ED",(0, (3, 5, 1, 5, 1, 5)),"firebrick" ,"MD",(0, (3, 1, 1, 1, 1, 1)),"royalblue" ],
            ["EQ",(0, (3, 5, 1, 5)),"maroon" ,"MQ",(0, (3, 1, 1, 1)), "navy"],
            ["EO",(0, (5, 1)),"indianred" ,"MO",(0, (5, 5)),"cornflowerblue"],
            ["E16",(0, (1, 1)),"lightcoral" ,"M16" ,(0, (1, 5)),"slateblue"]]

    if sizeParameter:
        var = "Size Parameter $ x=k a=\\frac{2 \pi n a}{\lambda}$"
    else:
        var = "$\lambda$"

    if multipoles is not None:
        n_num = len(multipoles)
        for n in range(1,n_num+1):
            if (n <= max_dis) and (n <= 4):
                ax.plot(waveLengths, multipoles[n-1][0], label = keys[n-1][0], linestyle = keys[n-1][1],color=keys[n-1][2])
                ax.plot(waveLengths, multipoles[n-1][1], label = keys[n-1][3], linestyle = keys[n-1][4],color=keys[n-1][5])
            elif (n <= max_dis):
                ax.plot(waveLengths, multipoles[n-1][0],color="red")
                ax.plot(waveLengths, multipoles[n-1][1],color="blue")

    ax.plot(waveLengths, scattering, c = "black", label = "Total",linewidth=0.8)
    ax.set_xlabel(var)
    if norm:
        ax.set_ylabel("$\sigma_{s}/a^2$ ")
    else:
        ax.set_ylabel("$\sigma_{s}$ ($m^{2}$)")
    ax.set_title("Scattering cross section of sphere of raduis, a = " + str(a) + prefix + " and $m=n_{1}/n=$" + str(m) + "\n calculated from Mie-Theroy displaying up to n = " + str(max_dis))
    #ax.set_title("Scattering cross section of a sphere showing the different scattering regimes")
    ax.legend()


if __name__ == "__main__":
    m = 3.55
    #m = np.sqrt(12.25)
    #perm = 12.25 + 1j*0.1
    #m = np.sqrt((np.abs(perm)+ perm.real)/2)
    #m = 3.9766
    
    
    a = 230
    MAX_dis = 2
    fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
    #waveLength, scattering, Multipoles = ScatteringCrossSection(m, a, 1000, 2000, multipole = True, norm = True, prefix="nm")
    waveLength, scattering, Multipoles = ScatteringCrossSection(m, a, 0.1, 10, multipole = True, norm = True, prefix="nm",sizeParameter = True)
    PlotScatteringCrossSection(ax, waveLength, scattering, multipoles = Multipoles,sizeParameter = True,max_dis = 4)




    #data = np.loadtxt('mieScatteringSi230nm_comsol_data_0.csv', delimiter=',')

    #comsol_waveLen = data[:,0]
    #comsol_scatter = data[:,1]



    #ax.plot(comsol_waveLen, comsol_scatter, label = "COMSOL",color="blue",linestyle="dashdot")
    #ax.legend()
    #plt.savefig('mieScattering_m{}_a{}nm.png'.format(m, round(a*10**9)))
    #plt.savefig('mieScattering_longboi.png')
    plt.show()