import tangos
import matplotlib.pyplot as plt

simname = "snapshots_200crit_h148mint"
timestep = tangos.get_timestep("snapshots_200crit_h148mint/h148.cosmo50PLK.6144g3HbwK1BH.004096")

m_star_halos = np.array([])
meanFeH_halos = np.array([])
meanOxH_halos = np.array([])
meanOxFe_halos = np.array([])

for halo in timestep.halos:
    m_star_halos = np.append(m_star_halos, halo['M_star'][0])
    meanFeH_halos = np.append(meanFeH_halos, halo['meanFeH'])
    meanOxH_halos = np.append(meanOxH_halos, halo['meanOxH'])
    meanOxFe_halos = np.append(meanOxFe_halos, halo['meanOxFe'])

fig, ax = plt.subplots()
ax.plot(m_star_halos, meanFeH_halos, 'ko')
ax.set_xscale('log')
ax.set_xlim(1e3, 1e10)
ax.set_ylim(-5.5,-0.5)
ax.set_xlabel("M_star")
ax.set_ylabel("[Fe/H]")
fig.savefig('/home/christenc/FeMstar.png')

fig, ax = plt.subplots()
ax.plot(m_star_halos, meanOxH_halos, 'ko')
ax.set_xscale('log')
ax.set_xlim(1e3, 1e10)
ax.set_ylim(-5.5,-0.5)
ax.set_xlabel("M_star")
ax.set_ylabel("[Ox/H]")
fig.savefig('/home/christenc/OxMstar.png')

fig, ax = plt.subplots()
ax.plot(m_star_halos, meanOxFe_halos, 'ko')
ax.set_xscale('log')
ax.set_xlim(1e3, 1e10)
ax.set_ylim(-0.25, 0.75)
ax.set_xlabel("M_star")
ax.set_ylabel("[Ox/Fe]")
fig.savefig('/home/christenc/OxFeMstar.png')

fig, ax = plt.subplots()
ax.plot(m_star_halos, meanOxH_halos, 'ko')
ax.set_xscale('log')
ax.set_xlim(1e3, 1e10)
ax.set_ylim(-5.5,-0.5)
ax.set_xlabel("M_star")
ax.set_ylabel("[Ox/H]")
fig.savefig('/home/christenc/OxMstar.png')

fig, ax = plt.subplots()
ax.plot(meanFeH_halos, meanOxFe_halos, 'ko')
ax.set_xlim(-5.5, -0.5)
ax.set_ylim(-0.25, 0.75)
ax.set_xlabel("[Fe/H]")
ax.set_ylabel("[Ox/Fe]")
fig.savefig('/home/christenc/OxFe_Fe.png')

