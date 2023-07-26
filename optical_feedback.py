# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 13:48:57 2023

@author: Malo Briend

Simulation of the coupled laser frequency versus the free-running 
laser frequency.

This code reproduces the figures of the article DOI: 10.1364/OE.431934.
You have to set the two boolean variables below : 
- 'non_res' defines if the non-resonant part of the equation is taken into 
account (default: True).
- 'experimental' defines if the experimental parameters are taken into account 
(default: False).
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

c = 299792458 # speed of light
rm2 = 0.99986 # external cavity mirror field reflection coefficient
r0 = 0.6 # exit facet power reflection coefficient
r02 = r0**2 # square of the exit facet power reflection coefficient
n = 3.5 # laser gain medium refractive index
a = 2*np.pi # switch to frequency
beta = 3e-5 # round-trip power attenuation
alpha = 2 # Henry-Factor
theta = np.arctan(alpha)# arctan of the Henry-Factor
Lc = 0.394 # Fabry-Pérot cavity length (m)
Ld = 0.001 # laser gain medium length (m)
F = 2*np.sqrt(rm2)/(1-rm2) # finesse
fsr = 3.9e8 # free spectral range
N = 100000 # number of points
freq_lim = 2e6 # frequency limit to create the array around one mode

non_res = True # add non-resonant part into the equation
experimental = False # takes experimental parameters instead of article's ones
if experimental == True:
    F = 30000
    beta = 1e-3
    Lc = 0.075
    rm2 = 0.9999
    fsr = 2e9
    
def coupling_cond(Lc=Lc, freq_lim = freq_lim):
    '''
    Find the integer to verify the coupling condition.
    Define the boolean variable limit to choose between a zoom or wide plot.
    Define the frequency range.
    '''
    tau_c = 2*Lc/c
    m = 80912
    wq = 2*m*np.pi/tau_c # wq around 1550 nm (in angular frequency)
    nu_q = wq/a
    tau_a = (2*m*np.pi - theta)/wq
    xmin = -2*np.pi*freq_lim + wq
    xmax = 2*np.pi*freq_lim + wq
    omega = np.linspace(xmin,xmax,N) # coupled angular frequency
    nu = omega/a
    return tau_c, tau_a, nu_q, omega, nu

def main(x, Lc=Lc, beta=beta, F=F):
    '''
    First, recover the coupling condition from the function coupling_cond().
    Calculate all prefactors and the part of the main equation.
    Deduce from that the locking range and Delta nu
    '''
    tau_c, tau_a, nu_q, omega, nu = coupling_cond(Lc, freq_lim)

    # feedback coupling rate pre-factor
    K = np.sqrt((1+alpha**2)*beta)*c*(1-r02)/(2*n*Ld*r0)
    K1 = K*np.sqrt(rm2)/(1-rm2) 
    K2 = K*np.sqrt(rm2)

    # non resonant term and resonant term
    nonres = K2*np.sin(x*tau_a+theta)
    res = K1*(np.sin(x*(tau_a+tau_c)+theta) - rm2*np.sin(x*tau_a+theta))\
        /(1 + F**2 * np.sin(x*tau_c/2)**2)
    omega_nonres = x - nonres
    nu_nonres = omega_nonres/a
    omega_res = x + res
    nu_res = omega_res/a

    if non_res == True:
        omega_coupled = x + res - nonres
    else:
        omega_coupled = x + res
    nu_coupled = omega_coupled/a
    max_lr = np.max(omega_res)
    min_lr = np.min(omega_res)
    locking_range = max_lr-min_lr

    D_omega_max = omega[np.where(omega_res == max_lr)[0]]
    D_omega_min = omega[np.where(omega_res == min_lr)[0]]
    Delta_omega = D_omega_max-D_omega_min
    Delta_nu = Delta_omega/a
    return nu_nonres, nu_coupled, nu_res, locking_range, Delta_nu

def free_running_freq(Lc=Lc):
    '''
    Reproduce the plot from the article.
    freq_lim = 400e6 to see 3 modes
    '''
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
    if experimental == False:
        freq_lim = 400e6
        ax2.set_xlim(-5e7,5e7)
        ax2.set_ylim(-2.2e6,2.2e6)
    else:
        freq_lim = 800e6
        ax2.set_xlim(-2e8,2e8)
        ax2.set_ylim(-2e6,2e6)
        
    nu_q, omega, nu = coupling_cond(freq_lim=freq_lim)[2:]
    nu_nonres, nu_coupled = main(omega, Lc=Lc)[:2]
    if non_res==True:
        ax1.plot(nu_nonres-nu_q, nu-nu_q, color='b', label='Only non resonant OF')
        ax2.plot(nu_nonres-nu_q, nu-nu_q, color='b')
    ax1.plot(nu-nu_q, nu-nu_q, linestyle='dashdot', color='k', label='Without OF')
    ax1.plot(nu_coupled-nu_q, nu-nu_q, linestyle='dashed', color='r', label='With total OF')
    ax2.plot(nu-nu_q, nu-nu_q, linestyle='dashdot', color='k')
    ax2.plot(nu_coupled-nu_q, nu-nu_q, linestyle='dashed', color='r')
    ax1.set_xlabel('Free running frequency (Hz)')
    ax2.set_xlabel('Free running frequency (Hz)')
    ax1.set_ylabel('Coupled frequency (Hz)')
    ax1.legend()
    plt.show()
    return fig

def beta_evolution(Lc=Lc):
    '''
    Compute the locking range for different beta values.
    '''
    if experimental == True:
        print('Error : please set the boolean variable experimental = False')
        return
    omega = coupling_cond(freq_lim=freq_lim)[3]
    list_beta = np.linspace(1e-7, 1e-3, 1000)
    locking_range1 = [main(omega, Lc=Lc, beta=i, F=10000)[3]*1e-9 for i in list_beta]
    locking_range2 = [main(omega, Lc=Lc, beta=i, F=20000)[3]*1e-9 for i in list_beta]
    locking_range3 = [main(omega, Lc=Lc, beta=i, F=100000)[3]*1e-9 for i in list_beta]

    fig, ax = plt.subplots()
    ax.set_xlabel(r'Feedback rate $\beta$')
    ax.set_ylabel('Locking range (GHz)')
    ax.plot(list_beta, locking_range1, color='darkorchid', label='F = 10000')
    ax.plot(list_beta, locking_range2, color='slateblue', label='F = 20000')
    ax.plot(list_beta, locking_range3,  color='cornflowerblue', label='F = 100000')
    ax.set_xscale('log')
    ax.legend()
    rect = Rectangle((5e-6, 0), 1e-4, 1.7, edgecolor='k', fc = 'lightgray')
    ax.add_patch(rect)
    ax.text(9.5e-6, 1.8, 'Usual rates')
    plt.show()
    return fig

def finesse_evolution(Lc=Lc):
    '''
    Compute the locking range for different finesse.
    Conclude on a limit if ratio > 1 = overlapping.
    '''
    if experimental == True:
        print('Error : please set the boolean variable experimental = False')
        return
    omega = coupling_cond(freq_lim=freq_lim)[3]
    list_finesse = np.linspace(3e3, 2e5, 1000)
    locking_range_fsr = [main(omega, Lc=Lc, F=i)[3]/fsr for i in list_finesse]
    Delta_nu = [main(omega, Lc=Lc, F=i)[4]*1e-3 for i in list_finesse]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
    ax1.set_xlabel('Finesse of the cavity')
    ax1.set_ylabel('Locking range/FSR')
    ax1.set_xscale('log')
    ax1.hlines(1, list_finesse[0], list_finesse[-1], linestyles='--', colors='k')
    ax1.plot(list_finesse, locking_range_fsr, color='slateblue')
    ax2.set_xlabel('Finesse of the cavity')
    ax2.set_ylabel(r'$\Delta \nu$ (kHz)')
    ax2.set_xscale('log')
    ax2.plot(list_finesse, Delta_nu, color='slateblue')
    plt.show()
    return fig
    
def upsample(x, y, ts_new=1):
    '''
    Add more point to known data.
    Variable ts_new is the new span between each point.
    '''
    fct = interpolate.interp1d(x, y, kind='cubic')
    x_upsamp = np.arange(x[0], x[-1], ts_new)
    y_upsamp = fct(x_upsamp)
    return x_upsamp, y_upsamp

def freq_jumps(Lc=Lc):
    '''
    freq_lim = 140e6 here to have the curve due to the non resonant term
    Define approximately the start and the end of the locking range. 
    On this part, add more point with the upsample() function.
    Concatenate the lists to have the new_data set.
    Find the two point of curvature where the jumps will begins.
    Create a constant line and find intersection between this line and free-running freq
    Get all the x and y coordinates and index.
    Concatenate the small list to one list and then feed the Airy function with that.
    Keep in mind that the axis are reversed.
    Hope it will works.
    EDIT : it works ._.
    '''
    if experimental == False:
        freq_lim = 140e6
    else:
        freq_lim = 800e6
    nu_q, omega, nu = coupling_cond(Lc, freq_lim=freq_lim)[2:]
    nu_nonres, nu_coupled = main(omega, Lc=Lc)[:2]

    list1 = nu-nu_q
    list2 = nu_coupled-nu_q
    start_plateau = 49970
    end_plateau = 51876

    list_ante_plateau_x = list1[:start_plateau]
    list_past_plateau_x = list1[end_plateau:N]
    list_ante_plateau_y = list2[:start_plateau]
    list_past_plateau_y = list2[end_plateau:N]
    x_upsamp, y_upsamp = upsample(list1[start_plateau:end_plateau], list2[start_plateau:end_plateau], ts_new=10)
    length_plateau = len(y_upsamp)

    # new data list with more point on the plateau of length 122041
    new_data_x = np.concatenate((list_ante_plateau_x, x_upsamp, list_past_plateau_x))
    new_data_y = np.concatenate((list_ante_plateau_y, y_upsamp, list_past_plateau_y))

    curve1 = np.max(new_data_y[:start_plateau])
    curve2 = np.max(new_data_y[start_plateau:start_plateau+length_plateau])
    i_min_jump2 = np.where(new_data_y == curve2)[0][0]

    line1 = [curve1]*len(new_data_x)
    line2 = [curve2]*len(new_data_x)
    idx1 = np.argwhere(np.diff(np.sign(line1-new_data_y))).flatten()
    idx2 = np.argwhere(np.diff(np.sign(line2-new_data_y))).flatten()
    
    y_min_jump1 = new_data_y[idx1][1]
    i_min_jump1 = np.where(new_data_y == y_min_jump1)[0][0]
    x_min_jump1 = new_data_x[i_min_jump1]
    y_max_jump1 = new_data_y[idx1+1][2]
    i_max_jump1 = np.where(new_data_y == y_max_jump1)[0][0]
    x_max_jump1 = new_data_x[i_max_jump1]

    y_min_jump2 = new_data_y[idx2][1]
    i_min_jump2 = np.where(new_data_y == y_min_jump2)[0][0]
    x_min_jump2 = new_data_x[i_min_jump2]
    y_max_jump2 = new_data_y[idx2+1][2]
    i_max_jump2 = np.where(new_data_y == y_max_jump2)[0][0]
    x_max_jump2 = new_data_x[i_max_jump2]

    line1_x = np.linspace(x_min_jump1, x_max_jump1, len(new_data_x))
    line2_x = np.linspace(x_min_jump2, x_max_jump2, len(new_data_x))   
    
    fig, ax = plt.subplots()
    list_jumps_y = np.concatenate((new_data_x[0:i_min_jump1], line1_x,\
                                    new_data_x[i_max_jump1:i_min_jump2],\
                                        line2_x, new_data_x[i_max_jump2:]))
    list_jumps_x = np.concatenate((new_data_y[0:i_min_jump1], line1,\
                                    new_data_y[i_max_jump1:i_min_jump2],\
                                        line2, new_data_y[i_max_jump2:]))

    new_data_x_MHz = [i*1e-6 for i in new_data_x]
    new_data_y_MHz = [i*1e-6 for i in new_data_y]
    list_jumps_x_MHz = [i*1e-6 for i in list_jumps_x]
    list_jumps_y_MHz = [i*1e-6 for i in list_jumps_y]
    
    ax.plot(new_data_y_MHz, new_data_x_MHz, label='Numerical simulation')
    ax.plot(list_jumps_x_MHz,list_jumps_y_MHz, color='r', label='Actual frequency jumps')
    ax.set_xlabel('Free running frequency (MHz)')
    ax.set_ylabel('Coupled frequency (MHz)')
    ax.legend()
    plt.show()
    return list_jumps_x, list_jumps_y

def transmission():
    '''
    Get the transmission of the cavity with the Airy function.
    '''
    def airy(x,R):
        phi = [(a*i*Lc)/c for i in x]
        return (1-R)**2 / ((1-R)**2 + 4*R*np.sin(phi)**2)

    list_jumps_x, list_jumps_y = freq_jumps()
    list_jumps_MHz = [i*1e-6 for i in list_jumps_x]
    airy_jumps = airy(list_jumps_y,rm2)

    nu_q, omega, nu = coupling_cond()[2:]
    new_xlist = [i-nu_q for i in nu]
    airy_free = airy(new_xlist,rm2)
    new_xlist_MHz = [i*1e-6 for i in new_xlist]

    fig, ax = plt.subplots()
    ax.plot(list_jumps_MHz, airy_jumps, color='r', label='With feedback')
    ax.plot(new_xlist_MHz, airy_free, color='b', label='Without feedback')
    ax.set_xlabel('Free running frequency (MHz)')
    ax.set_ylabel('Transmission')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    free_running_freq()
    transmission()
    beta_evolution()
    finesse_evolution()