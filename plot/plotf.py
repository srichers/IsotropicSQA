import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pylab as pylab
from multiprocessing import pool
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as colors
import h5py

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
orig_cmap = plt.get_cmap('gist_heat')
cmap = truncate_colormap(orig_cmap, 0, 0.85, 200)

clight = 3e10 # cm/s
h = 6.6260755e-27 # erg s
eV = 1.60218e-12 # erg

ns = 2
ng=25

mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2

E = (np.arange(0,25)+1)*4 #np.array([2.0, 2.43154, 2.95619, 3.59405, 4.36953, 5.31235, 6.45859, 7.85216, 9.54642, 11.6062, 14.1105, 17.1551, 20.8567, 25.3569, 30.8282, 37.48])

def line_color(ig,ng, cmap):
    return cmap(E[ig]/E[-1])

def column(a,b, state, group, filename):
    realindex = 1 + group*2*ns*ns*2 + state*ns*ns*2 + b*ns*2 + a*2
    imagindex = realindex+1
    realcol = np.genfromtxt(filename, usecols=(realindex))
    imagcol = np.genfromtxt(filename, usecols=(imagindex))[:len(realcol)]
    return realcol, imagcol

def isospin(f):
    fx = (f[:,:,0,1]+f[:,:,1,0])/2.
    fz = (f[:,:,0,0]-f[:,:,1,1])/2.
    ft = (f[:,:,0,0]+f[:,:,1,1])/2.
    return fx,fz,ft

def magnitude(r,i):
    return np.sqrt(r**2 + i**2)

def makeplot(location):
    filename = location+"/f.h5"
    f=h5py.File(filename,"r")
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8,12))
    t = np.array(f["r(cm)"]) / clight * 1e12
    fmatrixf = np.array(f["fmatrixf"])

    for ig in range(ng):
        r = fmatrixf[:,0,ig,0,1,0]
        axes[0].plot(t, fmatrixf[:,0,ig,0,0,0], '-', color=line_color(ig,ng,cmap), linewidth=1.5)
        p, = axes[1].plot(t, fmatrixf[:,0,ig,1,1,0], '-', color=line_color(ig,ng,cmap), linewidth=1.5)
        axes[2].plot(t, fmatrixf[:,0,ig,1,0,0]/fmatrixf[0,0,ig,1,0,0], '-', color=line_color(ig,ng,cmap), linewidth=1.5)
        axes[3].plot(t, fmatrixf[:,0,ig,1,0,1]/fmatrixf[0,0,ig,1,0,0], '-', color=line_color(ig,ng,cmap), linewidth=1.5)
        if(ig==0):
            p1=p
    for ig in range(ng):
        r = fmatrixf[:,1,ig,0,1,0]
        axes[0].plot(t, fmatrixf[:,1,ig,0,0,0], '--', color=line_color(ig,ng,cmap), linewidth=1.5)
        p, = axes[1].plot(t, fmatrixf[:,1,ig,1,1,0], '--', color=line_color(ig,ng,cmap), linewidth=1.5)
        axes[2].plot(t, fmatrixf[:,1,ig,1,0,0]/fmatrixf[0,1,ig,1,0,0], '--', color=line_color(ig,ng,cmap), linewidth=1.5)
        axes[3].plot(t, fmatrixf[:,1,ig,1,0,1]/fmatrixf[0,1,ig,1,0,0], '--', color=line_color(ig,ng,cmap), linewidth=1.5)
        if ig==0:
            p2=p
            axes[1].legend([p1,p2],[r"$f$",r"$\bar{f}$"],fontsize=22, ncol=2,loc='upper center',bbox_to_anchor=(0.5, 1), handletextpad=0, frameon=False)

    #axes[2].grid()
    axes[3].set_xlabel(r"$t\,(\mathrm{ps})$",fontsize=22)
    #axes[2].set_ylabel(r"$f_{e\mu}(t)/f_{e\mu}(0)$",fontsize=22)
    #axes[2].savefig("foffdiag.pdf",bbox_inches="tight")
    

    #####################################

        
    #if ig==0:
    #    axes[0].legend([p1,p2],[r"$f_{ee}$",r"$f_{\mu\mu}$"], fontsize=22,loc="upper right", ncol=2, frameon=False)
    #    axes[1].legend([p3,p4],[r"$\bar{f}_{ee}$",r"$\bar{f}_{\mu\mu}$"],fontsize=22,loc="upper right", ncol=2, frameon=False)

        
    # colorbar
    cax = fig.add_axes([.92, 0.1, 0.04, 0.8])
    norm = mpl.colors.Normalize(vmin=0, vmax=E[-1])
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,norm=norm,orientation='vertical')
    cb1.set_label(r'$h\nu\,(\mathrm{MeV})$')
    cax.minorticks_on()

    #axes[0].set_xlabel(r"$t\,\mathrm{(ms)}$",fontsize=22)
    #axes[0].set_ylabel(r"$f(t)/f(0)$",fontsize=22)
    ylabelx = -.15
    plt.text(ylabelx,0.5,r"$f_{ee}(t)$",horizontalalignment='center',verticalalignment='center', transform=axes[0].transAxes,rotation=90,fontsize=22)
    plt.text(ylabelx,0.5,r"$f_{\mu\mu}(t)$",horizontalalignment='center',verticalalignment='center', transform=axes[1].transAxes,rotation=90,fontsize=22)
    plt.text(ylabelx,0.5,r"$f_{(x)}(t)/f_{(x)}(0)$",horizontalalignment='center',verticalalignment='center', transform=axes[2].transAxes,rotation=90,fontsize=22)
    plt.text(ylabelx,0.5,r"$f_{(y)}(t)/f_{(x)}(0)$",horizontalalignment='center',verticalalignment='center', transform=axes[3].transAxes,rotation=90,fontsize=22)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.minorticks_on()
    axes[0].yaxis.set_major_locator(plt.MultipleLocator(.5))
    axes[1].yaxis.set_major_locator(plt.MultipleLocator(.2))
    axes[2].yaxis.set_major_locator(plt.MultipleLocator(.5))
    axes[3].yaxis.set_major_locator(plt.MultipleLocator(.5))
    #axes[0].set_ylim(0,1.001)
    #axes[1].set_ylim(0,1)
    axes[2].set_ylim(-1.4,1.4)
    axes[3].set_ylim(-1.4,1.4)
    for ax in axes:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        #ax.ticklabel_format(useOffset=False)
        ax.set_xlim(0,.68)
    
    plt.savefig(location+"/f_short.pdf",bbox_inches="tight")

    
for loc in ["."]:#["noosc_abs","noosc_fakepair","noosc_pair","noosc_fakebrems","noosc_escat","noosc_iscat","noosc_nointeract","noosc_nucscat"]:
    makeplot(loc)
