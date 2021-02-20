# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

# wg_EIM_profile - Effective Index Method - mode profile

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from matplotlib.patches import Rectangle

#thickness = 220e-9     # thickness of the silicon layer, [units: m]
#polarization = 'TE'  # polarization, 'TE' or 'TM'

lam, t, w, t_slab, n_core, n_clad, n_oxide, pts, M = 1.60e-6, 0.22e-6, 0.5e-6, 0, 3.47, 1.44, 1.44, 100, 2

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

    return x, TE_e, TE_h, TM_e, TM_h

    
#def draw_WG_vertical(M):
#	pP=get(gca,'Position');pPw=pP(3); 
#	pPc=pP(3)/2+pP(1); pP2=pPw/4/M;
#	annotation ('line',[pPc-pP2,pPc-pP2], [pP(2),pP(4)+pP(2)],'LineStyle','--');
#	annotation ('line',[pPc+pP2,pPc+pP2], [pP(2),pP(4)+pP(2)],'LineStyle','--');
#	axis tight; a=axis; axis ([a(1), a(2), a(3)*1.1, a(4)*1.1]);
#
#def draw_WG_horiz(M):
#	pP=get(gca,'Position');pPw=pP(4); 
#	pPc=pP(4)/2+pP(2); pP2=pPw/4/M;
#	annotation ('line',[pP(1),pP(3)+pP(1)], [pPc-pP2,pPc-pP2],'LineStyle','--');
#	annotation ('line',[pP(1),pP(3)+pP(1)], [pPc+pP2,pPc+pP2],'LineStyle','--');
#	axis tight; a=axis; axis ([a(1)*1.1, a(2)*1.1, a(3), a(4)]);

#==============================================================================
# Effective Index Method - mode profile
#==============================================================================
# find TE (TM) modes of slab waveguide (waveguide core and slab portions):
nTE, nTM, TEparam, TMparam = wg_1D_analytic(lam, t, n_oxide, n_core, n_clad)
if t_slab >= 0:
    nTE_slab, nTM_slab, TEparam_slab, TMparam_slab = wg_1D_analytic(lam, t_slab, n_oxide, n_core, n_clad)
else:
    nTE_slab = n_clad
    nTM_slab = n_clad

xslab, TE_Eslab, TE_Hslab, TM_Eslab, TM_Hslab = wg_1D_mode_profile(lam, t, n_oxide, n_core, n_clad, pts, M)
nx = np.hstack((nTE_slab[0]*np.ones(pts), nTE[0]*np.ones(pts), nTE_slab[0]*np.ones(pts)))

Fontsize=15

fig, ax1 = plt.subplots()
for i in np.arange(np.size(TE_Eslab[:,0])):
    ax1.plot(TE_Eslab[i]/np.max(np.max(TE_Eslab)), xslab*1e9, linewidth = 4)
ax1.set_ylabel('Height [nm]', fontsize = Fontsize)
ax1.set_xlabel('E-field (TE)', fontsize = Fontsize)
ax1.set_xticklabels([])
ax1.tick_params('x', colors='b')
# the second plot share y-axis with the first one
ax2 = ax1.twiny()
ax2.plot(nx, xslab*1e9, 'r--', linewidth = 0.5)
ax2.set_xlabel('Material Index', fontsize = Fontsize, color='r')
ax2.tick_params('x', colors='r')
ax2.set_xlim([np.min(nx)-0.05, np.max(nx)+0.05])
fig.tight_layout()
plt.show()

# TE-like modes of the etched waveguide (for fundamental slab mode)
#   solve for the "TM" modes:
nTE, nTM, TEparam, TMparam = wg_1D_analytic(lam, w, nTE_slab[0], nTE[0], nTE_slab[0])
neff_TEwg = nTM
xwg, TE_E_TEwg, TE_H_TEwg, TM_E_TEwg, TM_H_TEwg = wg_1D_mode_profile(lam, w, nTE_slab[0], nTE[0], nTE_slab[0], pts, M)
nx = np.hstack((nTE_slab[0]*np.ones(pts), nTE[0]*np.ones(pts), nTE_slab[0]*np.ones(pts)))

# Plot the data on with a left and right axes. Return the axes and line
fig, ax1 = plt.subplots()
for i in np.arange(np.size(TM_E_TEwg[:,0])):
    ax1.plot(xwg*1e9, TM_E_TEwg[i,:]/np.max(np.max(TM_E_TEwg)), linewidth = 4)
ax1.set_xlabel('Position [nm]', fontsize = Fontsize)
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('E-field (TM, TE-like mode)', fontsize = Fontsize, color='b')
ax1.tick_params('y', colors='b')
ax1.set_yticklabels([])
# the second plot share x-axis with the first one
ax2 = ax1.twinx()
ax2.plot(xwg*1e9, nx, 'r--')
ax2.set_ylabel('Slab Effective Index', fontsize = Fontsize, color='r')
ax2.tick_params('y', colors='r')
ax2.set_ylim([np.min(nx)-0.05, np.max(nx)+0.05])
fig.tight_layout()
plt.show()

# Plot the product of the two fields
plt.figure()
Exy = np.dot(TM_E_TEwg[0, :][:, np.newaxis], TE_Eslab[0, :][np.newaxis, :])
plt.contourf(xwg*1e9, xslab*1e9, np.abs(Exy.T)/np.max(np.max(Exy)).T, cmap=plt.cm.jet)
# colormap('jet')
plt.axis('equal')
plt.xlabel ('X (nm)', fontsize = Fontsize)
plt.ylabel ('Y (nm)', fontsize = Fontsize)
plt.title('Effective Index Method')
# Draw the waveguide:
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((-w/2*1e9, -t/2*1e9), w*1e9,t*1e9, alpha=1, facecolor='none', edgecolor='white', linewidth = 1))

if t_slab > 0:
    currentAxis.add_patch(Rectangle((-M*w*1e9, -t/2*1e9), (M-0.5)*w*1e9, t_slab*1e9, alpha=1, facecolor='none', edgecolor='white', linewidth = 1))
    currentAxis.add_patch(Rectangle((w/2*1e9, -t/2*1e9), (M-0.5)*1e9, t_slab*1e9, alpha=1, facecolor='none', edgecolor='white', linewidth = 1))
plt.show()

print ('neff_TEwg = ', neff_TEwg)
