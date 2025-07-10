# Charlotte Christensen

# 6/23/22
# Plot the rediduals SMHM relation for the marvel and Justice League runs, coloring points according to distance/tau_90

# SMHM for environment paper

#%run /home/christenc/Code/python/python_analysis/SMHM_v_distance_resid
import matplotlib as mpl
#mpl.use('tkagg') #Also can try mpl.use('Agg') #for using mpl over ssh    
#mpl.use('Agg')
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import pynbody
import math
import numpy as np
import socket
import matplotlib.gridspec as gridspec
import pandas as pd
import sys, os, glob, pickle
from scipy.interpolate import interp1d
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task
import Simpy
import tangos

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

outfile_base = '/home/christenc/Figures/marvel/marvelJL'
presentation = False
if presentation:
    outfile_base = outfile_base + '_pres'
    plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
    plt_width = 8 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 16
    dpi = 200
    markersize = 100
    ms_scale = 1
    lw = mpl.rcParams['lines.linewidth'] - 1
    edgewidth = 1
else:
    plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
    plt_width = 3.5 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 5
    dpi = 300
    markersize = 25
    ms_scale = 0.25
    lw = mpl.rcParams['lines.linewidth']
    edgewidth = 0.5

if (socket.gethostname() == "ozma.grinnell.edu"):
    dataprefix = '/home/christensen/Code/Datafiles/' 
else:
    dataprefix = '/home/christenc/Code/Datafiles/'

objs_pd = pd.read_csv(dataprefix + 'SMHM_env.csv')
objs_pd = objs_pd.set_index('halo_label')
objs_pd_resid = pd.read_csv(dataprefix + 'SMHM_env_resid.csv')
objs_pd_resid = objs_pd_resid.set_index('halo_label')
del objs_pd_resid['Mpeak']
del objs_pd_resid['Mstar_z0']
del objs_pd_resid['max_tindex']
del objs_pd_resid['min_dist']
del objs_pd_resid['type']
del objs_pd_resid['concentration']
objs_pd_comb = pd.concat([objs_pd,objs_pd_resid], join="inner", axis=1)

objs_pd_comb['p_size'] = np.ones(len(objs_pd_comb))*markersize
mask = objs_pd_comb[(objs_pd_comb['sim']=='h148') | (objs_pd_comb['sim']=='h229') | (objs_pd_comb['sim']=='h242') | (objs_pd_comb['sim']=='h329')].index
objs_pd_comb.loc[mask,'p_size'] = markersize*0.5

objs_pd_comb['max_tindex'] = objs_pd_comb['max_tindex'] #**(1/2)
objs_pd_comb['Mpeak_snap'] = objs_pd_comb['Mpeak_snap']*13.8/4096

panel_keys = ['massiveDist', 'min_dist', 'max_tindex', 'Mpeak_snap', 'tau90_vir', 'concentration']
panel_keys = ['massiveDist', 'min_dist', 'max_tindex', 'Mpeak_snap', 'tau90_vir', 'concent_highz']
#panel_names = ["D$_{\mathrm{massive}}$ [kpc]","Min(D$_{\mathrm{massive}}$) [kpc]","Max(Tidal Index)$^(1/2)$","Time of Peak M$_{\mathrm{vir}}$ [Gyr]",r"$\tau_{90}$ [Gyr]","Concentration at M$_{\mathrm{vir, peak}}$"]
panel_names = ["D$_{\mathrm{massive}}$ [kpc]","Min(D$_{\mathrm{massive}}$) [kpc]","Max(Tidal Index$)^(1/2)$","Time of Peak M$_{\mathrm{vir}}$ [Gyr]",r"$\tau_{90}$ [Gyr]","Concentration at high $z$"]
plts_log = [True, True, True, False, False, True]
n_panels = len(panel_names)

plt.close('all')
fig3 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
#fig3 = plt.figure(figsize=(plt_width,plt_width*2))
cNorm  = colors.LogNorm(vmin=1e8, vmax = 2e12) # Peak Virial mass
cmx = plt.get_cmap("winter_r") 
gs = fig3.add_gridspec(3)
axs = []
for count, (panel_key, panel_name, plt_log) in enumerate(zip(panel_keys[:3], panel_names[:3], plts_log[:3])):
    ax = fig3.add_subplot(gs[count])
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Central')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] == 13.8)])
    ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite') & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")
    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none", cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], alpha = 0.5)
    ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] < 13.8)]*1.5, marker = "*", alpha = 0.5)    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite') & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] < 13.8)]*1.5, marker = "*", alpha = 0.5)
        
    plt.axhline(0, linestyle = '--', color = 'k')
    min_x = min(objs_pd_comb[panel_key][objs_pd_comb[panel_key] > 0])
    x_fit = np.linspace(min_x, max(objs_pd_comb[panel_key]), 3)    
    if plt_log:
        ax.set_xscale('log')
        z = np.polyfit(np.log10(objs_pd_comb[panel_key])[(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['Mpeak_snap'] == 13.8)],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(np.log10(x_fit)), color = 'k')
        print(p(x_fit))        
    else:
        z = np.polyfit(objs_pd_comb[panel_key][objs_pd_comb['Mpeak_snap'] == 13.8], objs_pd_comb['resid1'][objs_pd_comb['Mpeak_snap'] == 13.8],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(x_fit), color = 'k')
        print(p(x_fit))
    ax.axis([min_x, max(objs_pd_comb[panel_key]), -5, 4])
    ax.set_xlabel(panel_name)
    ax.set_ylabel("Residual")
    ax.tick_params(axis = 'x', direction = 'in', top = False)
    ax.tick_params(axis = 'x', which = 'minor', direction = 'in', top = False)
    axs.append(ax)
    
fig3.tight_layout()
fig3.show()
sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
cbar = plt.colorbar(sm, location="right", ax = axs, aspect = 60) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig3.savefig(outfile_base + '_resid_partial.png',dpi = dpi,bbox_inches='tight')

panel_keys = ['massiveDist', 'Mpeak_snap', 'min_dist', 'tau90_vir', 'max_tindex', 'concentration']
panel_keys = ['massiveDist', 'Mpeak_snap', 'min_dist', 'tau90_vir', 'max_tindex', 'concent_highz']
panel_names = ["D$_{\mathrm{massive}}$ [kpc]","Time of Peak M$_{\mathrm{vir}}$ [Gyr]","Min(D$_{\mathrm{massive}}$) [kpc]",r"$\tau_{90}$ [Gyr]","Max(Tidal Index)","Concentration at M$_{\mathrm{vir, peak}}$"]
panel_names = ["D$_{\mathrm{massive}}$ [kpc]","Time of Peak M$_{\mathrm{vir}}$ [Gyr]","Min(D$_{\mathrm{massive}}$) [kpc]",r"$\tau_{90}$ [Gyr]","Max(Tidal Index)","Concentration at high $z$"]
plts_log = [True, False, True, False, True, True]
n_panels = len(panel_names)

plt.close('all')
#fig1 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*n_panels))
fig1 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11) # Peak Virial mass
cmx = plt.get_cmap("winter_r") 
gs = fig1.add_gridspec(ncols = 2, nrows = int(n_panels/2))
axs = []
for count, (panel_key, panel_name, plt_log) in enumerate(zip(panel_keys, panel_names, plts_log)):
    print(count, np.floor(count/2),count%2)
    if (count%2)==0:
        ax = fig1.add_subplot(gs[int(np.floor(count/2)),count%2])
    else:
        ax = fig1.add_subplot(gs[int(np.floor(count/2)),count%2], sharey = axs[-1])
    ax.scatter(objs_pd_comb[panel_key][objs_pd_comb['type']=='Central'], objs_pd_comb['resid1'][objs_pd_comb['type']=='Central'], c = objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'])
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*")    
    ax.scatter(objs_pd_comb[panel_key][objs_pd_comb['type']=='Satellite'], objs_pd_comb['resid1'][objs_pd_comb['type']=='Satellite'], c = objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = "*")
    plt.axhline(0, linestyle = '--', color = 'k')
    min_x = min(objs_pd_comb[panel_key][objs_pd_comb[panel_key] > 0])
    x_fit = np.linspace(min_x, max(objs_pd_comb[panel_key]), 3)    
    if plt_log:
        ax.set_xscale('log')
        z = np.polyfit(np.log10(objs_pd_comb[panel_key])[np.isfinite(np.log10(objs_pd_comb[panel_key]))], objs_pd_comb['resid1'][np.isfinite(np.log10(objs_pd_comb[panel_key]))],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(np.log10(x_fit)), color = 'k')
        slope = p[1]
        slope = (p(np.log10(x_fit[-1])) - p(np.log10(x_fit[0]))) #/(np.log10(x_fit[-1]) - np.log10(x_fit[0]))
        print(np.log10(x_fit[0]), np.log10(x_fit[-1]) ,p(np.log10(x_fit[0])),  p(np.log10(x_fit[-1])))
        print(p(x_fit))        
    else:
        z = np.polyfit(objs_pd_comb[panel_key], objs_pd_comb['resid1'],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(x_fit), color = 'k')
        slope = p[1]
        slope = (p(x_fit[-1]) - p(x_fit[0])) #/(x_fit[-1] - x_fit[0])
        print(x_fit[0], x_fit[-1],p(x_fit[0]),  p(x_fit[-1]))
        print(p(x_fit))
    if p[1] < 0:
        ax.text(0.06, 0.06, "{0:.3f}".format(slope), transform=ax.transAxes, fontsize = mpl.rcParams['axes.titlesize'])
    else:
        ax.text(0.8, 0.06, "{0:.3f}".format(slope),transform=ax.transAxes, fontsize = mpl.rcParams['axes.titlesize'])    
    ax.axis([min_x, max(objs_pd_comb[panel_key]), -5, 4])
    ax.set_xlabel(panel_name)
    if count%2 == 0:
        ax.set_ylabel("Residual")
    ax.tick_params(axis = 'x', direction = 'in', top = False)
    ax.tick_params(axis = 'x', which = 'minor', direction = 'in', top = False)
    axs.append(ax)
fig1.tight_layout()
fig1.show()
sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
cbar = plt.colorbar(sm, location="right", ax = axs, aspect = 60) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig1.savefig(outfile_base + '_resid2.png',dpi = dpi,bbox_inches='tight')

plt.close('all')
fig2 = plt.figure(15,figsize=(plt_width,plt_width*aspect_ratio))
gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
ax = fig2.add_subplot(gs[0])
ax.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'], c = objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'])
ax.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*")    
ax.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'], c = objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = "*")
ax.set_xscale('log')
ax.set_xlabel("Min(D$_{\mathrm{massive}}$) [kpc]")
ax.set_ylabel("Time of Peak M$_{\mathrm{vir}}$ [Gyr]")
fig2.tight_layout()
fig2.show()
cbar = plt.colorbar(sm, location="right", ax = [ax], aspect = 20) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig2.savefig(outfile_base + '_mindist_tmpeak.png',dpi = dpi,bbox_inches='tight')

plt.close('all')
fig3 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
#fig3 = plt.figure(figsize=(plt_width,plt_width*2))
cNorm  = colors.LogNorm(vmin=1e8, vmax = 2e12) # Peak Virial mass
cmx = plt.get_cmap("winter_r") 
gs = fig3.add_gridspec(3)
axs = []
for count, (panel_key, panel_name, plt_log) in enumerate(zip(panel_keys[:3], panel_names[:3], plts_log[:3])):
    ax = fig3.add_subplot(gs[count])
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Central')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] == 13.8)])
    ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite') & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")
    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none", cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], alpha = 0.5)
    ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] < 13.8)]*1.5, marker = "*", alpha = 0.5)    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite') & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] < 13.8)]*1.5, marker = "*", alpha = 0.5)
        
    plt.axhline(0, linestyle = '--', color = 'k')
    min_x = min(objs_pd_comb[panel_key][objs_pd_comb[panel_key] > 0])
    x_fit = np.linspace(min_x, max(objs_pd_comb[panel_key]), 3)    
    if plt_log:
        ax.set_xscale('log')
        z = np.polyfit(np.log10(objs_pd_comb[panel_key])[(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['Mpeak_snap'] == 13.8)],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(np.log10(x_fit)), color = 'k')
        print(p(x_fit))        
    else:
        z = np.polyfit(objs_pd_comb[panel_key][objs_pd_comb['Mpeak_snap'] == 13.8], objs_pd_comb['resid1'][objs_pd_comb['Mpeak_snap'] == 13.8],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(x_fit), color = 'k')
        print(p(x_fit))
    ax.axis([min_x, max(objs_pd_comb[panel_key]), -5, 4])
    ax.set_xlabel(panel_name)
    ax.set_ylabel("Residual")
    ax.tick_params(axis = 'x', direction = 'in', top = False)
    ax.tick_params(axis = 'x', which = 'minor', direction = 'in', top = False)
    axs.append(ax)
fig3.tight_layout()
fig3.show()
sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
cbar = plt.colorbar(sm, location="right", ax = axs, aspect = 60) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig3.savefig(outfile_base + '_resid_partial.png',dpi = dpi,bbox_inches='tight')


plt.close('all')
fig4 = plt.figure(figsize=(plt_width,plt_width*2))
#fig4 = plt.figure(figsize=(plt_width,plt_width*2))
cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11) # Peak Virial mass
cmx = plt.get_cmap("winter_r") 
gs = fig4.add_gridspec(n_panels)
axs = []
for count, (panel_key, panel_name, plt_log) in enumerate(zip(panel_keys, panel_names, plts_log)):
    ax = fig4.add_subplot(gs[count])
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central')], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central')], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Central')], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central')])
    #ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")    
    #ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite') & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")
    
    #ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none", cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], alpha = 0.5)
    ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5, marker = "*", alpha = 0.5)    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite')], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')]*1.5, marker = "*", alpha = 0.5)
        
    plt.axhline(0, linestyle = '--', color = 'k')
    min_x = min(objs_pd_comb[panel_key][objs_pd_comb[panel_key] > 0])
    x_fit = np.linspace(min_x, max(objs_pd_comb[panel_key]), 3)    
    if plt_log:
        ax.set_xscale('log')
        z = np.polyfit(np.log10(objs_pd_comb[panel_key])[(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['type'] == 'Central')], objs_pd_comb['resid1'][(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['type'] == 'Central')],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(np.log10(x_fit)), color = 'k')

        print(p(x_fit))        
    else:
        z = np.polyfit(objs_pd_comb[panel_key], objs_pd_comb['resid1'],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(x_fit), color = 'k')
        slope = p[1]
        slope = (p(x_fit[-1]) - p(x_fit[0])) #/(x_fit[-1] - x_fit[0])
        print(x_fit[0], x_fit[-1],p(x_fit[0]),  p(x_fit[-1]))
        print(p(x_fit))
    if p[1] < 0:
        ax.text(0.06, 0.06, "{0:.3f}".format(slope), transform=ax.transAxes, fontsize = mpl.rcParams['axes.titlesize'])
    else:
        ax.text(0.8, 0.06, "{0:.3f}".format(slope),transform=ax.transAxes, fontsize = mpl.rcParams['axes.titlesize'])    
    ax.axis([min_x, max(objs_pd_comb[panel_key]), -5, 4])
    ax.set_xlabel(panel_name)
    ax.set_ylabel("Residual")
    ax.tick_params(axis = 'x', direction = 'in', top = False)
    ax.tick_params(axis = 'x', which = 'minor', direction = 'in', top = False)
    axs.append(ax)
sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
fig4.tight_layout()
fig4.show()
cbar = plt.colorbar(sm, location="right", ax = axs, aspect = 60) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig4.savefig(outfile_base + '_resid.png',dpi = dpi,bbox_inches='tight')

plt.close('all')
fig2 = plt.figure(15,figsize=(plt_width,plt_width*aspect_ratio))
gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
ax = fig2.add_subplot(gs[0])
ax.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'], c = objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'])
ax.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*")    
ax.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'], c = objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = "*")
ax.set_xscale('log')
ax.set_xlabel("Min(D$_{\mathrm{massive}}$) [kpc]")
ax.set_ylabel("Time of Peak M$_{\mathrm{vir}}$ [Gyr]")
fig2.tight_layout()
fig2.show()
cbar = plt.colorbar(sm, location="right", ax = [ax], aspect = 20) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig2.savefig(outfile_base + '_mindist_tmpeak.png',dpi = dpi,bbox_inches='tight')


plt.close('all')
fig4 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
#fig4 = plt.figure(figsize=(plt_width,plt_width*2))
cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11) # Peak Virial mass
cmx = plt.get_cmap("winter_r") 
gs = fig4.add_gridspec(ncols = 2, nrows = int(n_panels/2))
axs = []
for count, (panel_key, panel_name, plt_log) in enumerate(zip(panel_keys, panel_names, plts_log)):
    print(count, np.floor(count/2),count%2)
    if (count%2)==0:
        ax = fig4.add_subplot(gs[int(np.floor(count/2)),count%2])
    else:
        ax = fig4.add_subplot(gs[int(np.floor(count/2)),count%2], sharey = axs[-1])    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central')], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central')], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Central')], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central')])
    #ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")    
    #ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite') & (objs_pd_comb['Mpeak_snap'] == 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], c = objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)], cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')  & (objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, marker = "*")
    
    #ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], facecolor = "none", c = "none", cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & (objs_pd_comb['Mpeak_snap'] < 13.8)], alpha = 0.5)
    ax.scatter(objs_pd_comb[panel_key][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))], objs_pd_comb['resid1'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5, marker = "*", alpha = 0.5)    
    ax.scatter(objs_pd_comb[panel_key][(objs_pd_comb['type']=='Satellite')], objs_pd_comb['resid1'][(objs_pd_comb['type']=='Satellite')], facecolor = "none", c = "none",  cmap = cmx, norm = cNorm, edgecolor = 'k', linewidths = edgewidth, s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite')]*1.5, marker = "*", alpha = 0.5)
        
    plt.axhline(0, linestyle = '--', color = 'k')
    min_x = min(objs_pd_comb[panel_key][objs_pd_comb[panel_key] > 0])
    x_fit = np.linspace(min_x, max(objs_pd_comb[panel_key]), 3)    
    if plt_log:
        ax.set_xscale('log')
        z = np.polyfit(np.log10(objs_pd_comb[panel_key])[(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['type'] == 'Central')], objs_pd_comb['resid1'][(np.isfinite(np.log10(objs_pd_comb[panel_key]))) & (objs_pd_comb['type'] == 'Central')],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(np.log10(x_fit)), color = 'k')
        slope = p[1]
        slope = (p(np.log10(x_fit[-1])) - p(np.log10(x_fit[0]))) #/(np.log10(x_fit[-1]) - np.log10(x_fit[0]))
        print(np.log10(x_fit[0]), np.log10(x_fit[-1]) ,p(np.log10(x_fit[0])),  p(np.log10(x_fit[-1])))        
        #print(p(x_fit))        
    else:
        z = np.polyfit(objs_pd_comb[panel_key][objs_pd_comb['type'] == 'Central'], objs_pd_comb['resid1'][objs_pd_comb['type'] == 'Central'],1)
        p = np.poly1d(z)
        ax.plot(x_fit, p(x_fit), color = 'k')
        slope = p[1]
        slope = (p(x_fit[-1]) - p(x_fit[0])) #/(x_fit[-1] - x_fit[0])
        print(x_fit[0], x_fit[-1],p(x_fit[0]),  p(x_fit[-1]))
        #print(p(x_fit))
    print(p, p(0))
    if p[1] < 0:
        ax.text(0.06, 0.06, "{0:.3f}".format(slope), transform=ax.transAxes, fontsize = mpl.rcParams['axes.titlesize'])
    else:
        ax.text(0.8, 0.06, "{0:.3f}".format(slope),transform=ax.transAxes, fontsize = mpl.rcParams['axes.titlesize'])       
    ax.axis([min_x, max(objs_pd_comb[panel_key]), -5, 4])
    ax.set_xlabel(panel_name)
    if count%2 == 0:
        ax.set_ylabel("Residual")
    ax.tick_params(axis = 'x', direction = 'in', top = False)
    ax.tick_params(axis = 'x', which = 'minor', direction = 'in', top = False)
    axs.append(ax)
fig4.tight_layout()
fig4.show()
sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
cbar = plt.colorbar(sm, location="right", ax = axs, aspect = 60) #pad=0.04)
cbar.set_label(r'M$_{\mathrm{vir, peak}} $ [M$_\odot$]') #, rotation=0)
fig4.savefig(outfile_base + '_resid_partial2.png',dpi = dpi,bbox_inches='tight')
