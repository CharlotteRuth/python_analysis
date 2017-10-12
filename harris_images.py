s = pynbody.load('h239.cosmo50cmb.3072g14HMbwK.00512')
outfilebase = 'h239.cosmo50cmb.3072g14HMbwK.00512'
s.physical_units()
h = s.halos()

fig = plt.figure(0)
pynbody.analysis.angmom.faceon(h[1])
sph.image(h[1].d,qty="rho",units="Msol pc^-2",width=500,cmap="Greys",vmin = 1e-1,vmax= 1e3,qtytitle='',title="Dark Matter")
fig.savefig(outfilebase + '_dm_cb.png')

fig = plt.figure(0)
sph.image(h[1].g,qty="rho",units="Msol pc^-2",width=500,cmap="gist_heat_r",vmin = 1e-1,vmax= 1e2,qtytitle='',title="Gas")
fig.savefig(outfilebase + '_gas_cb.png')

fig = plt.figure(0)
sph.image(h[1].s,qty="rho",units="Msol pc^-2",width=500,cmap="Blues",vmin = 1e-6,vmax= 1e4,qtytitle='',title="Stars")
fig.savefig(outfilebase + '_stars_cb.png')


fig = plt.figure(0)
pynbody.analysis.angmom.faceon(h[1])
sph.image(h[1].d,qty="rho",units="Msol pc^-2",width=500,cmap="Greys",vmin = 1e-1,vmax= 1e3,qtytitle='',show_cbar=False,title="Dark Matter")
fig.savefig(outfilebase + '_dm.png')

fig = plt.figure(0)
sph.image(h[1].g,qty="rho",units="Msol pc^-2",width=500,cmap="gist_heat_r",vmin = 1e-1,vmax= 1e2,qtytitle='',show_cbar=False,title="Gas")
fig.savefig(outfilebase + '_gas.png')

fig = plt.figure(0)
sph.image(h[1].s,qty="rho",units="Msol pc^-2",width=500,cmap="Blues",vmin = 1e-6,vmax= 1e4,qtytitle='',show_cbar=False,title="Stars")
fig.savefig(outfilebase + '_stars.png')
