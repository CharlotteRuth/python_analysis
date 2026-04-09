import tangos


tangos.core.init_db('/home/christenc/Databases/Marvel_BN_N10.db') # rom25_dwarf_zooms.db
tangos.all_simulations() # list all simulations
simname = "snapshots_200crit_h148mint" # Stored on emu; Define the simulation you want to look at
simname = "snapshots_200crit_h148" # Stored on quirm
tangos.get_simulation(simname).timesteps # list timesteps
tangos.get_simulation(simname).timesteps[-1].halos[0:100] # list first 100 halos in final timestep

tangos.get_simulation(simname).timesteps[-1].halos[0].keys() # list available keys

for halo in tangos.get_simulation(simname).timesteps[-1].halos:
    if halo['M_star'] > 0: #halo['M_star'][0] > 0:
        print(halo, halo['Mvir'], halo['M_star'])

