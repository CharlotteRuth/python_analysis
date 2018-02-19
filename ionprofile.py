def make_rs(im):
    xsize, ysize = np.shape(im)
    x = np.arange(-xsize/2, xsize/2)
    y = np.arange(-ysize/2, ysize/2)
    xs, ys = np.meshgrid(x,y)
    # 2.0 because using 2 kpc pixels
    return 2.0*np.sqrt(xs**2 + ys**2)


s = pynbody.load(tfile)
h = s.halos()

halo = h.load_copy(int(halo_num))
halo.physical_units()
pynbody.analysis.angmom.faceon(halo)
oviif = pynbody.analysis.ionfrac.calculate(halo.gas,ion='ovi')
halo.gas['oviden'] = halo.gas['rho']*halo.gas['OxMassFrac']*oviif
halo.gas['oviden'].units = halo.gas['rho'].units


oviim = pynbody.plot.image(halo.gas,qty='oviden', clear=False,
                   units='16 m_p cm^-2', width=1000,
                   vmin=12,vmax=17, log=True)
set_xlabel('x [kpc]')
rs = make_rs(oviim)

pynbody.plot.hist2d(np.log10(rs.flatten() + 1e-6),np.log10(oviim.flatten()),xrange=[0.3,3],
                    yrange=[10,17], cmap=plt.cm.gist_yarg)
plt.ylim(9,14)
plt.xlim(0.3,3)
plt.xlabel(r'log$_{10}$($\rho$ [kpc])')
plt.show()