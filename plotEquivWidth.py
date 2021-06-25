
#Run with
#%run /home/christensen/Code/python/python_analysis/plotEquivWidth.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/plotEquivWidth.py
#ipython --pylab

import numpy as np
import pynbody
import sys, os, glob, pickle
import socket
import pandas as pd
import matplotlib.pylab as plt

def pickle_read(file):

    objs = []
    f=open(file, 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs.append(p)
            #objs.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()

    return pd.DataFrame(objs)

if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        outprefix = '/home/christenc/Figures/marvel/'
        dataprefix = '/home/christenc/Code/Datafiles/'        
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        dataprefix = '/home/christensen/Code/Datafiles/'
   
    presentation = True
    if presentation:
        outbase = outprefix + 'marvel_pres_'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 12
        labelsize = 18
        dpi = 100
    else:
        outbase = outprefix + 'marvel'
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        labelsize = 5
        dpi = 300       

    #Cpt Marvel
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_cm = []
    f=open(tfile + '.data', 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs_cm.append(p)
        except EOFError:
            break
        
    f.close()
    objs_pd_cm = pd.DataFrame(objs_cm)

    #Elektra
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_e = []
    f=open(tfile + '.data', 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs_e.append(p)
        except EOFError:
            break
    f.close()
    objs_pd_e = pd.DataFrame(objs_e)

    #Rogue
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_r = []
    f=open(tfile + '.data', 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs_r.append(p)
        except EOFError:
            break
        
    f.close()
    objs_pd_r = pd.DataFrame(objs_r)

    #Storm
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_s = []
    f=open(tfile + '.data', 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs_s.append(p)
        except EOFError:
            break
        
    f.close()
    objs_pd_s = pd.DataFrame(objs_s)

    cm_1_label = 'CM1' #'{:.1e}'.format(objs_pd_cm['mstar'][0]) + r' M$_\odot$, CM 1' #'Cpt. Marvel, 1'
    e_1_label = 'E1' #'{:.1e}'.format(objs_pd_e['mstar'][0]) + r' M$_\odot$, E1' #'Elektra, 1'
    e_2_label = 'E2'#{:.1e}'.format(objs_pd_e['mstar'][1]) + r' M$_\odot$, E2' #'Elektra, 2'
    r_1_label = 'R1' #'{:.1e}'.format(objs_pd_s['mstar'][0]) + r' M$_\odot$, R1' #'Rogue, 1'
    r_3_label = 'R3' #{:.1e}'.format(objs_pd_s['mstar'][1]) + r' M$_\odot$, R3' #'Rogue, 3'
    s_1_label = 'S1' #{:.1e}'.format(objs_pd_r['mstar'][0]) + r' M$_\odot$, S1' #'Storm, 1'
    s_2_label = 'S2' #{:.1e}'.format(objs_pd_r['mstar'][1]) + r' M$_\odot$, S2' #'Storm, 2'
    zorder_cm_1 = 7
    zorder_e_1 = 3
    zorder_e_2 = 6
    zorder_r_1 = 2
    zorder_r_3 = 4
    zorder_s_1 = 1
    zorder_s_2 = 5
    markers = ['d','d','d','d','d','d','d']    
    colormap = 'tab20'
    obscolor = 'k'
    obsmark = 'o'
    alpha = 0.5
    colors_cm_1_f = 'none'
    colors_r_1_f = 'none'
    colors_e_1_f = 'none'
    colors_s_1_f = 'none'
    colors_r_3_f = 'none'
    colors_e_2_f = 'none'
    colors_s_2_f = 'none'
    colors_cm_1_e = 'darkviolet'
    colors_e_1_e = 'orange'
    colors_e_2_e = 'b'
    colors_r_1_e = 'r'
    colors_r_3_e = 'g' #'m'
    colors_s_1_e = 'deeppink' #'purple'
    colors_s_2_e = 'c'
 
#Find all line of sight files in the directory
    path = prefix + "elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "elektra.cosmo25cmb.4096g5HbwK1BH.004096"
    s = pynbody.load(path + tfile)
    h_dummy = s.halos(dummy = True)
    rvir_e_1 = h_dummy[1].properties['Rvir']

    objs_pd_e_1 =  pickle_read(path + tfile+'_1_CGM.data')
    objs_pd_e_1_v =  pickle_read(path + tfile+'_1_CGMvoight.data')
    
    rvir_e_2 = h_dummy[2].properties['Rvir']
    objs_pd_e_2 =  pickle_read(path + tfile+'_2_CGM.data')
    objs_pd_e_2_v =  pickle_read(path + tfile+'_2_CGMvoight.data')
    
    path = prefix + "rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "rogue.cosmo25cmb.4096g5HbwK1BH.004096"
    s = pynbody.load(path + tfile)
    h_dummy = s.halos(dummy = True)
    rvir_r_1 = h_dummy[1].properties['Rvir']
    objs_pd_r_1 =  pickle_read(path + tfile+'_1_CGM.data')
    objs_pd_r_1_v =  pickle_read(path + tfile+'_1_CGMvoight.data')

    rvir_r_3 = h_dummy[3].properties['Rvir']
    objs_pd_r_3 =  pickle_read(path + tfile+'_3_CGM.data')
    objs_pd_r_3_v =  pickle_read(path + tfile+'_3_CGMvoight.data')

    path = prefix + "storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "storm.cosmo25cmb.4096g5HbwK1BH.004096"
    s = pynbody.load(path + tfile)
    h_dummy = s.halos(dummy = True)
    rvir_s_1 = h_dummy[1].properties['Rvir']
    objs_pd_s_1 =  pickle_read(path + tfile+'_1_CGM.data')
    objs_pd_s_1_v =  pickle_read(path + tfile+'_1_CGMvoight.data')

    rvir_s_2 = h_dummy[2].properties['Rvir']
    objs_pd_s_2 =  pickle_read(path + tfile+'_2_CGM.data')
    objs_pd_s_2_v =  pickle_read(path + tfile+'_2_CGMvoight.data')

    path = prefix + "cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096"
    cm = pynbody.load(path + tfile)
    h_dummy = cm.halos(dummy = True)
    rvir_cm_1 = h_dummy[1].properties['Rvir']
    objs_pd_cm_1 =  pickle_read(path + tfile+'_1_CGM.data')
    objs_pd_cm_1_v =  pickle_read(path + tfile+'_1_CGMvoight.data')

    f = open(dataprefix+'Bordoloi2014.txt', 'r')
    bordoloidata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 19 and columns[0] != 'QSOName':
            source = {}
            source['qsoname'] = columns[0]
            source['galaxy'] = columns[1]
            source['alpha'] = columns[2]
            source['delta'] = columns[3]
            source['zsys'] = float(columns[4])
            source['L'] = float(columns[5])
            source['logMstar'] = float(columns[6])
            source['R'] = float(columns[7])
            source['Rvir'] = float(columns[8])
            source['logsSFR'] = float(columns[9])
            source['logsSFr_limit'] = columns[10]
            source['logN'] = float(columns[11])
            source['logN_err'] = float(columns[12])
            source['logN_limit'] = columns[13]
            source['Wr'] = float(columns[14])
            source['Wr_err'] = float(columns[15])
            source['Wr_limit'] = columns[16]
            source['phi'] = float(columns[17])
            source['CIR3sigma'] = float(columns[18])
            bordoloidata.append(source)
    f.close()
    bordoloidata = pd.DataFrame(bordoloidata)
    exact_wr = bordoloidata['Wr_limit'] == "no"
    upper_wr = bordoloidata['Wr_limit'] == "upper"

    f = open(dataprefix+'Johnson2015.txt', 'r')
    johnsondata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 26 and columns[0] != 'ID':
            source = {}
            source['ID'] = columns[0]
            source['RA'] = columns[1]
            source['dec'] = columns[2]
            source['z_gal'] = columns[3]
            source['Mr'] = float(columns[4])
            source['logMstar'] = float(columns[5])
            source['d'] = float(columns[6])
            source['Rhalo'] = float(columns[7])
            source['Wr_HI'] = float(columns[8])
            source['Wr_HI_err'] = float(columns[9])
            source['Wr_HI_limit'] = columns[10]
            source['Wr_SiII'] = float(columns[11])
            source['Wr_SiII_err'] = float(columns[12])
            source['Wr_SiII_limit'] = columns[13]
            source['Wr_SiIII'] = float(columns[14])
            source['Wr_SiIII_err'] = float(columns[15])
            source['Wr_SiIII_limit'] = columns[16]            
            source['Wr_SiIV'] = float(columns[17])
            source['Wr_SiIV_err'] = float(columns[18])
            source['Wr_SiIV_limit'] = columns[19]
            source['Wr_CIV'] = float(columns[20])
            source['Wr_CIV_err'] = float(columns[21])
            source['Wr_CIV_limit'] = columns[22]            
            source['Wr_OVI'] = float(columns[23])
            source['Wr_OVI_err'] = float(columns[24])
            source['Wr_OVI_limit'] = columns[25]
            johnsondata.append(source)
    f.close()
    johnsondata = pd.DataFrame(johnsondata)

    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[johnsondata["Wr_HI_limit"] == "no"]['d']/johnsondata[johnsondata['Wr_HI_limit'] == "no"]['Rhalo'],johnsondata[johnsondata['Wr_HI_limit'] == "no"]['Wr_HI'],johnsondata[johnsondata['Wr_HI_limit'] == "no"]['Wr_HI_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    johnson_plt = plt.errorbar(johnsondata[johnsondata["Wr_HI_limit"] == "upper"]['d']/johnsondata[johnsondata['Wr_HI_limit'] == "upper"]['Rhalo'],johnsondata[johnsondata['Wr_HI_limit'] == "upper"]['Wr_HI'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_HI_limit'] == "upper"]['Wr_HI'],fmt=None,barsabove = True,color = obscolor,ecolor=obscolor,uplims=True,zorder = 10,label = '_nolegend_') #,capsize=3,mew=0)        
    cm_1_plt = plt.scatter(objs_pd_cm_1['b']/rvir_cm_1,objs_pd_cm_1['eq_width_hi'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label)
    plt.scatter(objs_pd_cm_1_v['b']/rvir_cm_1,objs_pd_cm_1_v['eq_width_hi'],marker = 'o',facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label)
    e_1_plt = plt.scatter(objs_pd_e_1['b']/rvir_e_1,objs_pd_e_1['eq_width_hi'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_1_v['b']/rvir_e_1,objs_pd_e_1_v['eq_width_hi'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    e_2_plt = plt.scatter(objs_pd_e_2['b']/rvir_e_2,objs_pd_e_2['eq_width_hi'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_e_2_v['b']/rvir_e_2,objs_pd_e_2_v['eq_width_hi'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    r_1_plt = plt.scatter(objs_pd_r_1['b']/rvir_r_1,objs_pd_r_1['eq_width_hi'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_1_v['b']/rvir_r_1,objs_pd_r_1_v['eq_width_hi'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    r_3_plt = plt.scatter(objs_pd_r_3['b']/rvir_r_3,objs_pd_r_3['eq_width_hi'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_r_3_v['b']/rvir_r_3,objs_pd_r_3_v['eq_width_hi'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    s_1_plt = plt.scatter(objs_pd_s_1['b']/rvir_s_1,objs_pd_s_1['eq_width_hi'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_1_v['b']/rvir_s_1,objs_pd_s_1_v['eq_width_hi'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    s_2_plt = plt.scatter(objs_pd_s_2['b']/rvir_s_2,objs_pd_s_2['eq_width_hi'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_s_2_v['b']/rvir_s_2,objs_pd_s_2_v['eq_width_hi'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(1,1,r'H I Ly$\alpha$',size = labelsize,color = 'k')
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([0.1, 3, 0.001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.HI_Rvir.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[johnsondata["Wr_HI_limit"] == "no"]['d'],johnsondata[johnsondata['Wr_HI_limit'] == "no"]['Wr_HI'],johnsondata[johnsondata['Wr_HI_limit'] == "no"]['Wr_HI_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_HI_limit"] == "upper"]['d'],johnsondata[johnsondata['Wr_HI_limit'] == "upper"]['Wr_HI'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_HI_limit'] == "upper"]['Wr_HI'],fmt=None,barsabove = True,color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')    
    plt.scatter(objs_pd_cm_1['b'],objs_pd_cm_1['eq_width_hi'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b'],objs_pd_e_1['eq_width_hi'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b'],objs_pd_e_2['eq_width_hi'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b'],objs_pd_r_1['eq_width_hi'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b'],objs_pd_r_3['eq_width_hi'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b'],objs_pd_s_1['eq_width_hi'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b'],objs_pd_s_2['eq_width_hi'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b'],objs_pd_cm_1_v['eq_width_hi'],marker = 'o',facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b'],objs_pd_e_1_v['eq_width_hi'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b'],objs_pd_e_2_v['eq_width_hi'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b'],objs_pd_r_1_v['eq_width_hi'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b'],objs_pd_r_3_v['eq_width_hi'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b'],objs_pd_s_1_v['eq_width_hi'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b'],objs_pd_s_2_v['eq_width_hi'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho$ [kpc]')
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(100,1,r'H I Ly$\alpha$',size = labelsize,color = 'k')
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([8, 300, 0.001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.HI_kpc.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['d']/johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['Rhalo'],johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['Wr_CIV'],johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['Wr_CIV_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_CIV_limit"] == "upper"]['d']/johnsondata[johnsondata['Wr_CIV_limit'] == "upper"]['Rhalo'],johnsondata[johnsondata['Wr_CIV_limit'] == "upper"]['Wr_CIV'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_CIV_limit'] == "upper"]['Wr_CIV'],fmt="_",color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')
    bordoloi_plt = plt.errorbar(bordoloidata[exact_wr]['R']/bordoloidata[exact_wr]['Rvir'],bordoloidata[exact_wr]['Wr']/1000,yerr=bordoloidata[exact_wr]['Wr_err']/1000,color = "grey",fmt = obsmark,zorder = 10, label = 'Bordoloi+ 2014')
    plt.errorbar(bordoloidata[upper_wr]['R']/bordoloidata[upper_wr]['Rvir'],bordoloidata[upper_wr]['Wr']/1000,xerr=None,yerr=.2*bordoloidata[upper_wr]['Wr']/1000,fmt="_",color = "grey",ecolor="grey",uplims=True,label = '_nolegend_') #,capsize=3,mew=0)    
    plt.scatter(objs_pd_cm_1['b']/rvir_cm_1,objs_pd_cm_1['eq_width_civ'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b']/rvir_e_1,objs_pd_e_1['eq_width_civ'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b']/rvir_e_2,objs_pd_e_2['eq_width_civ'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b']/rvir_r_1,objs_pd_r_1['eq_width_civ'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b']/rvir_r_3,objs_pd_r_3['eq_width_civ'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b']/rvir_s_1,objs_pd_s_1['eq_width_civ'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b']/rvir_s_2,objs_pd_s_2['eq_width_civ'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b']/rvir_cm_1,objs_pd_cm_1_v['eq_width_civ'],marker = 'o',facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b']/rvir_e_1,objs_pd_e_1_v['eq_width_civ'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b']/rvir_e_2,objs_pd_e_2_v['eq_width_civ'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b']/rvir_r_1,objs_pd_r_1_v['eq_width_civ'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b']/rvir_r_3,objs_pd_r_3_v['eq_width_civ'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b']/rvir_s_1,objs_pd_s_1_v['eq_width_civ'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b']/rvir_s_2,objs_pd_s_2_v['eq_width_civ'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(1,1,r'C IV 1548',size = labelsize,color = 'k')  
    plt.yscale('log')
    plt.xscale('log')
    #plt.legend()
    plt.axis([0.1, 3, 0.001, 2])
    plt.savefig(outbase + '.CIV_Rvir.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['d'],johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['Wr_CIV'],johnsondata[(johnsondata["Wr_CIV_limit"] == "no") & (johnsondata["Wr_CIV"] != -1)]['Wr_CIV_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_CIV_limit"] == "upper"]['d'],johnsondata[johnsondata['Wr_CIV_limit'] == "upper"]['Wr_CIV'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_CIV_limit'] == "upper"]['Wr_CIV'],fmt="_",color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')
    bordoloi_plt = plt.errorbar(bordoloidata[exact_wr]['R'],bordoloidata[exact_wr]['Wr']/1000,yerr=bordoloidata[exact_wr]['Wr_err']/1000,color = "grey",fmt = obsmark,zorder = 10, label = 'Bordoloi+ 2014')
    plt.errorbar(bordoloidata[upper_wr]['R'],bordoloidata[upper_wr]['Wr']/1000,xerr=None,yerr=.2*bordoloidata[upper_wr]['Wr']/1000,fmt="_",color = "grey",ecolor="grey",uplims=True,label = '_nolegend_') #,capsize=3,mew=0)    
    plt.scatter(objs_pd_cm_1['b'],objs_pd_cm_1['eq_width_civ'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b'],objs_pd_e_1['eq_width_civ'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b'],objs_pd_e_2['eq_width_civ'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b'],objs_pd_r_1['eq_width_civ'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b'],objs_pd_r_3['eq_width_civ'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b'],objs_pd_s_1['eq_width_civ'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b'],objs_pd_s_2['eq_width_civ'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b'],objs_pd_cm_1_v['eq_width_civ'],marker = 'o',facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b'],objs_pd_e_1_v['eq_width_civ'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b'],objs_pd_e_2_v['eq_width_civ'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b'],objs_pd_r_1_v['eq_width_civ'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b'],objs_pd_r_3_v['eq_width_civ'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b'],objs_pd_s_1_v['eq_width_civ'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b'],objs_pd_s_2_v['eq_width_civ'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho$ [kpc]')
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(100,1,r'C IV 1548',size = labelsize,color = 'k')  
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([8, 300, 0.001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.CIV_kpc.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['d']/johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['Rhalo'],johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['Wr_SiIV'],johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['Wr_SiIV_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_SiIV_limit"] == "upper"]['d']/johnsondata[johnsondata['Wr_SiIV_limit'] == "upper"]['Rhalo'],johnsondata[johnsondata['Wr_SiIV_limit'] == "upper"]['Wr_SiIV'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_SiIV_limit'] == "upper"]['Wr_SiIV'],fmt="_",color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')    
    plt.scatter(objs_pd_cm_1['b']/rvir_cm_1,objs_pd_cm_1['eq_width_siiv'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b']/rvir_e_1,objs_pd_e_1['eq_width_siiv'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b']/rvir_e_2,objs_pd_e_2['eq_width_siiv'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b']/rvir_r_1,objs_pd_r_1['eq_width_siiv'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b']/rvir_r_3,objs_pd_r_3['eq_width_siiv'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b']/rvir_s_1,objs_pd_s_1['eq_width_siiv'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b']/rvir_s_2,objs_pd_s_2['eq_width_siiv'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b']/rvir_cm_1,objs_pd_cm_1_v['eq_width_siiv'],marker = "o",facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b']/rvir_e_1,objs_pd_e_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b']/rvir_e_2,objs_pd_e_2_v['eq_width_siiv'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b']/rvir_r_1,objs_pd_r_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b']/rvir_r_3,objs_pd_r_3_v['eq_width_siiv'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b']/rvir_s_1,objs_pd_s_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b']/rvir_s_2,objs_pd_s_2_v['eq_width_siiv'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(1,1,r'Si IV 1393',size = labelsize,color = 'k')
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([0.1, 3, 0.0001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.SiIV_Rvir.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['d'],johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['Wr_SiIV'],johnsondata[(johnsondata["Wr_SiIV_limit"] == "no") & (johnsondata["Wr_SiIV"] != -1)]['Wr_SiIV_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_SiIV_limit"] == "upper"]['d'],johnsondata[johnsondata['Wr_SiIV_limit'] == "upper"]['Wr_SiIV'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_SiIV_limit'] == "upper"]['Wr_SiIV'],fmt="_",color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')    
    plt.scatter(objs_pd_cm_1['b'],objs_pd_cm_1['eq_width_siiv'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b'],objs_pd_e_1['eq_width_siiv'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b'],objs_pd_e_2['eq_width_siiv'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b'],objs_pd_r_1['eq_width_siiv'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b'],objs_pd_r_3['eq_width_siiv'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b'],objs_pd_s_1['eq_width_siiv'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b'],objs_pd_s_2['eq_width_siiv'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b'],objs_pd_cm_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b'],objs_pd_e_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b'],objs_pd_e_2_v['eq_width_siiv'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b'],objs_pd_r_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b'],objs_pd_r_3_v['eq_width_siiv'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b'],objs_pd_s_1_v['eq_width_siiv'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b'],objs_pd_s_2_v['eq_width_siiv'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho$ [kpc]')    
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(100,1,r'Si IV 1393',size = labelsize,color = 'k')   
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([8, 300, 0.0001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.SiIV_kpc.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['d']/johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['Rhalo'],johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['Wr_OVI'],johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['Wr_OVI_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_OVI_limit"] == "upper"]['d']/johnsondata[johnsondata['Wr_OVI_limit'] == "upper"]['Rhalo'],johnsondata[johnsondata['Wr_OVI_limit'] == "upper"]['Wr_OVI'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_OVI_limit'] == "upper"]['Wr_OVI'],fmt="_",color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')
    plt.scatter(objs_pd_cm_1['b']/rvir_cm_1,objs_pd_cm_1['eq_width_ovi'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b']/rvir_e_1,objs_pd_e_1['eq_width_ovi'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b']/rvir_e_2,objs_pd_e_2['eq_width_ovi'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b']/rvir_r_1,objs_pd_r_1['eq_width_ovi'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b']/rvir_r_3,objs_pd_r_3['eq_width_ovi'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b']/rvir_s_1,objs_pd_s_1['eq_width_ovi'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b']/rvir_s_2,objs_pd_s_2['eq_width_ovi'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b']/rvir_cm_1,objs_pd_cm_1_v['eq_width_ovi'],marker = "o",facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b']/rvir_e_1,objs_pd_e_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b']/rvir_e_2,objs_pd_e_2_v['eq_width_ovi'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b']/rvir_r_1,objs_pd_r_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b']/rvir_r_3,objs_pd_r_3_v['eq_width_ovi'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b']/rvir_s_1,objs_pd_s_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b']/rvir_s_2,objs_pd_s_2_v['eq_width_ovi'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'$W_r$ [$\AA$]')
    plt.text(1,1,r'O VI 1031',size = labelsize,color = 'k')   
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([0.1, 3, 0.001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.OVI_Rvir.png',dpi = dpi)
    
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    johnson_plt = plt.errorbar(johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['d'],johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['Wr_OVI'],johnsondata[(johnsondata["Wr_OVI_limit"] == "no") & (johnsondata["Wr_OVI"] != -1)]['Wr_OVI_err'],color = obscolor,fmt = obsmark,zorder = 10, label = 'Johnson+ 2017')
    plt.errorbar(johnsondata[johnsondata["Wr_OVI_limit"] == "upper"]['d'],johnsondata[johnsondata['Wr_OVI_limit'] == "upper"]['Wr_OVI'],xerr=None,yerr=0.2*johnsondata[johnsondata['Wr_OVI_limit'] == "upper"]['Wr_OVI'],fmt="_",color = obscolor,ecolor=obscolor,uplims=True,label = '_nolegend_')    
    plt.scatter(objs_pd_cm_1['b'],objs_pd_cm_1['eq_width_ovi'],marker = markers[0],facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1['b'],objs_pd_e_1['eq_width_ovi'],marker = markers[1],facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2['b'],objs_pd_e_2['eq_width_ovi'],marker = markers[2],facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1['b'],objs_pd_r_1['eq_width_ovi'],marker = markers[3],facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3['b'],objs_pd_r_3['eq_width_ovi'],marker = markers[4],facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1['b'],objs_pd_s_1['eq_width_ovi'],marker = markers[5],facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2['b'],objs_pd_s_2['eq_width_ovi'],marker = markers[6],facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)
    plt.scatter(objs_pd_cm_1_v['b'],objs_pd_cm_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_cm_1_f, edgecolors = colors_cm_1_e, alpha = alpha,label = cm_1_label, zorder = zorder_cm_1)
    plt.scatter(objs_pd_e_1_v['b'],objs_pd_e_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_e_1_f, edgecolors = colors_e_1_e, alpha = alpha,label = e_1_label, zorder = zorder_e_1)
    plt.scatter(objs_pd_e_2_v['b'],objs_pd_e_2_v['eq_width_ovi'],marker = 'o',facecolors = colors_e_2_f, edgecolors = colors_e_2_e, alpha = alpha,label = e_2_label, zorder = zorder_e_2)
    plt.scatter(objs_pd_r_1_v['b'],objs_pd_r_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_r_1_f, edgecolors = colors_r_1_e, alpha = alpha,label = r_1_label, zorder = zorder_r_1)
    plt.scatter(objs_pd_r_3_v['b'],objs_pd_r_3_v['eq_width_ovi'],marker = 'o',facecolors = colors_r_3_f, edgecolors = colors_r_3_e, alpha = alpha,label = r_3_label, zorder = zorder_r_3)
    plt.scatter(objs_pd_s_1_v['b'],objs_pd_s_1_v['eq_width_ovi'],marker = 'o',facecolors = colors_s_1_f, edgecolors = colors_s_1_e, alpha = alpha,label = s_1_label, zorder = zorder_s_1)
    plt.scatter(objs_pd_s_2_v['b'],objs_pd_s_2_v['eq_width_ovi'],marker = 'o',facecolors = colors_s_2_f, edgecolors = colors_s_2_e, alpha = alpha,label = s_2_label, zorder = zorder_s_2)    
    plt.xlabel(r'$\rho$ [kpc]')    
    plt.ylabel(r'$W_r$ [$\AA$]')   
    plt.text(100,1,r'O VI 1031',size = labelsize,color = 'k')       
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.axis([8, 300, 0.001, 2])
#    plt.legend(fontsize = legendsize)
    plt.savefig(outbase + '.OVI_kpc.png',dpi = dpi)
