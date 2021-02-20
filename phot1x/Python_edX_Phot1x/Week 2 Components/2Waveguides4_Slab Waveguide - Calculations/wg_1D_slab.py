# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

#==============================================================================
# Python Script to solve for the modes of the slab waveguide,
# calculate the effective indices
# plot the E-field profile for all the modes
# the user provides:
# - thickness of the silicon layer
# - desired polarization to plot
#==============================================================================

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from matplotlib.patches import Rectangle

thickness = 220e-9     # thickness of the silicon layer, [units: m]
polarization = 'TE'  # polarization, 'TE' or 'TM'
lam, t, n1, n2, n3, pts, M = 1.55e-6, thickness, 1.444, 3.47, 1.444, 100, 3

def TE_eq(b0,k0,n1,n2,n3,t):
    h0 = np.sqrt( (n2*k0)**2 - b0**2 )
    q0 = np.sqrt( b0**2 - (n1*k0)**2 )
    p0 = np.sqrt( b0**2 - (n3*k0)**2 )
    # the objective is to find zeroes of te0 and tm0
    te0 = np.tan( h0*t ) - (p0+q0)/h0/(1-p0*q0/h0**2)
    return te0, h0, q0, p0

def TM_eq(b0,k0,n1,n2,n3,t):
    h0 = np.sqrt( (n2*k0)**2 - b0**2 )
    q0 = np.sqrt( b0**2 - (n1*k0)**2 )
    p0 = np.sqrt( b0**2 - (n3*k0)**2 )
    pbar0 = (n2/n3)**2*p0
    qbar0 = (n2/n1)**2*q0
    tm0 = np.tan( h0*t ) - h0*(pbar0+qbar0)/(h0**2-pbar0*qbar0)
    return tm0, h0, q0, p0
    
#==============================================================================
# wg_1D_analytic.m - Analytic solution of a slab waveguide
# by Lumerical Solutions, http://www.lumerical.com/mode_online_help/slab_wg.m
# modified by Lukas Chrostowski, 2012
# See Yariv Photonics book, Chapter 3
# finds the TE and TM effective indices of a 3-layer slab waveguide
# 
# usage:
#  - get effective indices for supported modes:
#  [nTE, nTM] = wg_1D_analytic (1.55e-6, 0.22e-6, 1.444, 3.47, 1.444)
#  - optional, for plotting: TEparam,TMparam: h, q, p parameters of the mode.
#==============================================================================

def wg_1D_analytic(lam, t, n1, n2, n3):
    k0 = 2*np.pi/lam
    b0 = np.linspace(np.max((n1, n3))*k0, n2*k0, 1000)  # k0*n3 less than b less than k0*n2
    b0 = b0[:len(b0)-1]
    te0 = TE_eq(b0,k0,n1,n2,n3,t)[0]
    tm0 = TM_eq(b0,k0,n1,n2,n3,t)[0]
    
    # TE
    intervals = (te0>=0).astype(int)-(te0<0).astype(int)
    izeros = np.nonzero(np.diff(intervals) < 0)[0]
    X0=np.vstack((b0[izeros], b0[izeros+1])).T
    nzeros, scrap=X0.shape
    
    nTE = np.zeros(nzeros)
    TEparam = np.zeros((nzeros, 4))
    for i in np.arange(nzeros):
        def vTE_eq(x):
            return TE_eq(x,k0,n1,n2,n3,t)
        
        nTE[i] = fsolve(lambda x: TE_eq(x,k0,n1,n2,n3,t)[0], X0[i,0]) / k0
        TEparam[i,0],TEparam[i,1],TEparam[i,2],TEparam[i,3] = TE_eq(nTE[i]*k0,k0,n1,n2,n3,t)

    nTE = nTE[::-1]
    TEparam = TEparam[::-1, :]
    
    # TM
    intervals = (tm0>=0).astype(int)-(tm0<0).astype(int)
    izeros = np.nonzero(np.diff(intervals) < 0)[0]
    X0=np.vstack((b0[izeros], b0[izeros+1])).T
    nzeros, scrap=X0.shape
    
    nTM = np.zeros(nzeros)
    TMparam = np.zeros((nzeros, 4))
    for i in np.arange(nzeros):
        def vTM_eq(x):
            return TM_eq(x,k0,n1,n2,n3,t)
        
        nTM[i] = fsolve(lambda x: TM_eq(x,k0,n1,n2,n3,t)[0], X0[i,0]) / k0
        TMparam[i,0],TMparam[i,1],TMparam[i,2],TMparam[i,3] = TM_eq(nTM[i]*k0,k0,n1,n2,n3,t)

    if nzeros >= 0:
        nTM = nTM[::-1]
        TMparam = TMparam[::-1, :]
    else:
        nTM = []
    
    return nTE,nTM,TEparam,TMparam
   
#==============================================================================
# wg_1D_mode_profile.m - Calculate the 1D mode profile of a slab waveguide
# by Lukas Chrostowski, 2012
# See Yariv Photonics book, Chapter 3.2
# - function returns mode profiles for TE and TM modes (E, H components)
# usage, e.g.:
#  [x, TE_e, TE_h, TM_e, TM_h, nTE, nTM] = wg_1D_mode_profile (1.55e-6, 0.22e-6, 1.444, 3.47, 1.444, 100, 4)
#  plot (x, TE_e) 
#  plot (x, TM_e)
#==============================================================================

def wg_1D_mode_profile(lam, t, n1, n2, n3, pts, M):
    nTE,nTM,TEparam,TMparam = wg_1D_analytic(lam, t, n1, n2, n3)
    x1 = np.linspace(-M*t, -t/2, pts)
    x2 = np.linspace(-t/2, t/2, pts)
    x3 = np.linspace(t/2, M*t, pts)
    x = np.hstack((x1, x2, x3))
    nx = np.hstack((n1*np.ones(pts), n2*np.ones(pts), n3*np.ones(pts)))
    
    # constants
    mu0 = 4*np.pi*1e-7
    epsilon0 = 8.85e-12
    eta = np.sqrt(mu0/epsilon0)
    c = 3e8

    TE_e = np.zeros([np.size(nTE), np.size(nx)])
    TE_h = np.zeros([np.size(nTE), np.size(nx)])
    TM_e = np.zeros([np.size(nTE), np.size(nx)])
    TM_h = np.zeros([np.size(nTE), np.size(nx)])
    for i in np.arange(np.size(nTE)):
        h = TEparam[i,1]
        q = TEparam[i,2]
        p = TEparam[i,3]
        beta = 2*np.pi*nTE[i]/lam
        C = 2*h*np.sqrt(2*np.pi*c/lam*mu0 / (beta*(t+1/q+1/p)*(h**2+q**2))) # normalize to 1W
        # n1, n2, n3 regions
        TE_e[i,:] = C * np.hstack(( np.exp(q*(x1+t/2)), (np.cos(h*(x2+t/2)) + q/h*np.sin(h*(x2+t/2))), (np.cos(h*t) + q/h*np.sin(h*t))*np.exp(-p*(x3-t/2)) ))
    TE_h = TE_e*nx/eta

    for i in np.arange(np.size(nTM)):
        h = TMparam[i,1]
        q = TMparam[i,2]
        p = TMparam[i,3]
        qb = n2**2/n1**2*q
#        pb = n2**2/n3**2*p
        beta = 2*np.pi*nTM[i]/lam
        temp = (qb**2+h**2)/qb**2*(t/n2**2+(q**2+h**2)/(qb**2+h**2)/n1**2/q + (p**2+h**2)/(p**2+h**2)/n3**2/p)
        C = 2*np.sqrt( 2*np.pi*c/lam*epsilon0 / (beta*temp) ) # normalize to 1W
        TM_h[i,:] = C * np.hstack(( h/qb*np.exp(q*(x1+t/2)), (h/qb*np.cos(h*(x2+t/2))+np.sin(h*(x2+t/2))), (h/qb*np.cos(h*t)+np.sin(h*t))*np.exp(-p*(x3-t/2)) ))
    TM_e=TM_h/nx*eta

    return x, TE_e, TE_h, TM_e, TM_h, nTE, nTM

#def wg_1D_slab():
x, TE_e, TE_h, TM_e, TM_h, nTE, nTM = wg_1D_mode_profile(1.55e-6, thickness, 1.444, 3.47, 1.444, 100, 3)
plt.figure()
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((-thickness*1e9/2, -1), thickness*1e9, 2, alpha=1, facecolor=[.8, .8, .8]))
# g = rectangle('Position',[,],'FaceColor',[.8 .8 .8])
if polarization == 'TE':  
    for i in np.arange(np.size(TE_e[:,0])):
        plt.plot(x*1e9, TE_e[i,:]/np.max(np.max(TE_e)), linewidth = 3)
    plt.title('TE polarized mode(s)')
else:
    for i in np.arange(np.size(TM_e[:,0])):
        plt.plot (x*1e9, TM_e[i,:]/np.max(np.max(TM_e)), linewidth = 3)
    plt.title('TM polarized mode(s)')

plt.xlabel('Position - perpendicular to wafer [nm]')
plt.ylabel('E-Field Amplitude')
# ax.XTick = [floor(min(x)*1e9/100)*100:200:floor(max(x)*1e9/100)*100];
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.grid()
plt.box()
plt.show()

print ("Effective index value(s) of the TE mode(s):  ", nTE)
print ("Effective index value(s) of the TM mode(s):  ", nTM)


    

